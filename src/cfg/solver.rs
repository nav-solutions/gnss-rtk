//! Solver configuration preset

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Base (zenith referenced) pseudo-range standard deviation, in meters.
const DEFAULT_CODE_SIGMA_M: f64 = 1.0;

/// Base (zenith referenced) phase-range standard deviation, in meters.
const DEFAULT_PHASE_SIGMA_M: f64 = 3.0E-3;

/// Lower bound applied to sin(elevation) when forming the elevation weighting,
/// to avoid the singularity at the horizon (caps the model around 2.9°).
const MIN_SIN_ELEV: f64 = 0.05;

const fn default_max_gdop() -> f64 {
    5.0
}

const fn default_postfit_denoising() -> f64 {
    1000.0
}

const fn default_open_loop() -> bool {
    false
}

const fn default_weighting() -> Weighting {
    Weighting::None
}

/// [Weighting] strategy used to build the measurement weight matrix `W`,
/// where each measurement contributes `1/σ²` on the diagonal.
#[derive(Default, Clone, Copy, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Weighting {
    /// Uniform weighting: every measurement is trusted equally (`W = I`).
    /// This reproduces an ordinary (unweighted) least-squares solution and is
    /// the current default, preserving legacy behavior.
    #[cfg_attr(feature = "serde", serde(alias = "none", alias = "None"))]
    #[default]
    None,

    /// Elevation dependent weighting: `σ²(E) = σ₀² / sin²(E)`.
    /// Low-elevation measurements are down-weighted, as they suffer larger
    /// tropospheric residuals, multipath and antenna gain roll-off. This
    /// mirrors the tropospheric mapping function and helps most when the
    /// measurement pool is of heterogeneous quality.
    #[cfg_attr(feature = "serde", serde(alias = "elevation", alias = "Elevation"))]
    Elevation,
}

#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct SolverOpts {
    /// GDOP threshold to invalidate ongoing GDOP
    #[cfg_attr(feature = "serde", serde(default = "default_max_gdop"))]
    pub max_gdop: f64,

    /// GNSS-RTK allows the user to open the navigation filter loop,
    /// which is totally forbidden in roaming applications.
    #[cfg_attr(feature = "serde", serde(default = "default_open_loop"))]
    pub open_loop: bool,

    /// Possible extra denoising filter, at the expense
    /// of more processing time. The configuration is the denoising factor.
    /// 1000 for x1000 improvement attempt.
    #[cfg_attr(feature = "serde", serde(default = "default_postfit_denoising"))]
    pub postfit_denoising: f64,

    /// Measurement [Weighting] model used to build the weight matrix.
    #[cfg_attr(feature = "serde", serde(default = "default_weighting"))]
    pub weighting: Weighting,
}

impl Default for SolverOpts {
    fn default() -> Self {
        Self {
            max_gdop: default_max_gdop(),
            open_loop: default_open_loop(),
            postfit_denoising: default_postfit_denoising(),
            weighting: default_weighting(),
        }
    }
}

impl SolverOpts {
    /// Parameter settings recommended for static ultra precise applications
    pub fn static_preset() -> Self {
        Self {
            max_gdop: 3.0,
            open_loop: default_open_loop(),
            postfit_denoising: default_postfit_denoising(),
            weighting: default_weighting(),
        }
    }

    /// Returns the pseudo-range standard deviation (m) at `elevation_deg`,
    /// under the configured [Weighting] model.
    pub(crate) fn code_sigma_m(&self, elevation_deg: Option<f64>) -> f64 {
        self.measurement_sigma_m(DEFAULT_CODE_SIGMA_M, elevation_deg)
    }

    /// Returns the phase-range standard deviation (m) at `elevation_deg`,
    /// under the configured [Weighting] model.
    pub(crate) fn phase_sigma_m(&self, elevation_deg: Option<f64>) -> f64 {
        self.measurement_sigma_m(DEFAULT_PHASE_SIGMA_M, elevation_deg)
    }

    /// Returns the measurement standard deviation (m) for an observation of
    /// base (zenith) deviation `base_sigma_m` seen at `elevation_deg`, under
    /// the configured [Weighting] model. [Weighting::None] yields a unit
    /// deviation (uniform weighting); when elevation is unknown under an
    /// elevation model, the base deviation is returned.
    fn measurement_sigma_m(&self, base_sigma_m: f64, elevation_deg: Option<f64>) -> f64 {
        match self.weighting {
            Weighting::None => 1.0,
            Weighting::Elevation => match elevation_deg {
                Some(elev_deg) => {
                    let sin_elev = elev_deg.to_radians().sin().max(MIN_SIN_ELEV);
                    base_sigma_m / sin_elev
                },
                None => base_sigma_m,
            },
        }
    }
}

#[cfg(test)]
mod test {
    use super::{DEFAULT_CODE_SIGMA_M, MIN_SIN_ELEV, SolverOpts, Weighting};

    #[test]
    fn elevation_weighting() {
        let opts = SolverOpts {
            weighting: Weighting::Elevation,
            ..Default::default()
        };

        // zenith: deviation equals the base value
        let zenith = opts.code_sigma_m(Some(90.0));
        assert!((zenith - DEFAULT_CODE_SIGMA_M).abs() < 1.0E-9);

        // 30°: sin(30°)=0.5 -> deviation doubles
        let e30 = opts.code_sigma_m(Some(30.0));
        assert!((e30 - DEFAULT_CODE_SIGMA_M / 0.5).abs() < 1.0E-9);

        // low elevation is down-weighted (larger sigma than at zenith)
        assert!(opts.code_sigma_m(Some(5.0)) > zenith);

        // phase is far more precise than code
        assert!(opts.phase_sigma_m(Some(90.0)) < opts.code_sigma_m(Some(90.0)));

        // unknown elevation falls back to the base deviation
        assert!((opts.code_sigma_m(None) - DEFAULT_CODE_SIGMA_M).abs() < 1.0E-9);
    }

    #[test]
    fn uniform_weighting() {
        let opts = SolverOpts {
            weighting: Weighting::None,
            ..Default::default()
        };

        // uniform weighting => unit deviation regardless of elevation
        assert_eq!(opts.code_sigma_m(Some(10.0)), 1.0);
        assert_eq!(opts.phase_sigma_m(Some(90.0)), 1.0);
        assert_eq!(opts.code_sigma_m(None), 1.0);
    }

    #[test]
    fn horizon_guard() {
        let opts = SolverOpts {
            weighting: Weighting::Elevation,
            ..Default::default()
        };

        // sigma stays finite at/below the horizon thanks to MIN_SIN_ELEV
        let sigma = opts.code_sigma_m(Some(0.0));
        assert!(sigma.is_finite());
        assert!((sigma - DEFAULT_CODE_SIGMA_M / MIN_SIN_ELEV).abs() < 1.0E-9);
    }
}
