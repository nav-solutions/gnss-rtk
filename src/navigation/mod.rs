pub mod solutions;
pub use solutions::{PVTSolution, PVTSolutionType};

mod filter;

pub use filter::{Filter, FilterState};

use log::debug;
use std::collections::HashMap;

use crate::{
    bias::{Bias, IonosphereBias, RuntimeParam as BiasRuntimeParams, TropoModel, TroposphereBias},
    candidate::Candidate,
    cfg::Config,
    prelude::{Error, Method, SV},
};

use nalgebra::{
    base::dimension::{U4, U8},
    OMatrix, OVector,
};
use nyx::cosmic::SPEED_OF_LIGHT;

/// SV Navigation information
#[derive(Debug, Clone, Default)]
pub struct SVInput {
    /// SV azimuth angle in degrees
    pub azimuth: f64,
    /// SV elevation angle in degrees
    pub elevation: f64,
    /// Ionospheric bias in meters of delay
    pub iono_bias: Bias,
    /// Tropospheric bias in meters of delay
    pub tropo_bias: Bias,
}

/// Navigation Input
#[derive(Debug, Clone)]
pub struct Input {
    /// Measurement vector
    pub y: OVector<f64, U8>,
    /// NAV Matrix
    pub g: OMatrix<f64, U8, U8>,
    /// Weight Diagonal Matrix
    pub w: OMatrix<f64, U8, U8>,
    /// SV dependent data
    pub sv: HashMap<SV, SVInput>,
}

/// Navigation Output
#[derive(Debug, Clone, Default)]
pub struct Output {
    /// Time Dilution of Precision
    pub tdop: f64,
    /// Geometric Dilution of Precision
    pub gdop: f64,
    /// Position Dilution of Precision
    pub pdop: f64,
    /// Q covariance matrix
    pub q: OMatrix<f64, U8, U8>,
    /// Filter state
    pub state: FilterState,
}

impl Output {
    pub(crate) fn q_covar4x4(&self) -> OMatrix<f64, U4, U4> {
        OMatrix::<f64, U4, U4>::new(
            self.q[(0, 0)],
            self.q[(0, 1)],
            self.q[(0, 2)],
            self.q[(0, 3)],
            self.q[(1, 0)],
            self.q[(1, 1)],
            self.q[(1, 2)],
            self.q[(1, 3)],
            self.q[(2, 0)],
            self.q[(2, 1)],
            self.q[(2, 2)],
            self.q[(2, 3)],
            self.q[(3, 0)],
            self.q[(3, 1)],
            self.q[(3, 2)],
            self.q[(3, 3)],
        )
    }
}

impl Input {
    /// Forms new Navigation Input
    pub fn new(
        apriori: (f64, f64, f64),
        apriori_geo: (f64, f64, f64),
        cfg: &Config,
        cd: &[Candidate],
        iono_bias: &IonosphereBias,
        tropo_bias: &TroposphereBias,
    ) -> Result<Self, Error> {
        let mut y = OVector::<f64, U8>::zeros();
        let mut g = OMatrix::<f64, U8, U8>::zeros();
        let mut sv = HashMap::<SV, SVInput>::with_capacity(cd.len());
        /*
         * Compensate for ARP (if possible)
         */
        let apriori = match cfg.arp_enu {
            Some(offset) => (
                apriori.0 + offset.0,
                apriori.1 + offset.1,
                apriori.2 + offset.2,
            ),
            None => apriori,
        };

        let (x0, y0, z0) = apriori;

        for i in 0..8 {
            let mut sv_input = SVInput::default();

            let index = if i >= cd.len() {
                if cfg.sol_type == PVTSolutionType::TimeOnly {
                    0
                } else {
                    i - cd.len()
                }
            } else {
                i
            };

            let state = cd[index].state.ok_or(Error::UnresolvedState)?;
            let clock_corr = cd[index].clock_corr.to_seconds();

            let (azimuth, elevation) = (state.azimuth, state.elevation);
            sv_input.azimuth = azimuth;
            sv_input.elevation = elevation;

            let (sv_x, sv_y, sv_z) = (state.position[0], state.position[1], state.position[2]);
            let rho = ((sv_x - x0).powi(2) + (sv_y - y0).powi(2) + (sv_z - z0).powi(2)).sqrt();
            let (x_i, y_i, z_i) = ((x0 - sv_x) / rho, (y0 - sv_y) / rho, (z0 - sv_z) / rho);

            g[(i, 0)] = x_i;
            g[(i, 1)] = y_i;
            g[(i, 2)] = z_i;
            g[(i, 3)] = 1.0_f64;

            let mut models = 0.0_f64;

            if cfg.modeling.sv_clock_bias {
                models -= clock_corr * SPEED_OF_LIGHT;
            }
            if let Some(delay) = cfg.externalref_delay {
                models -= delay * SPEED_OF_LIGHT;
            }

            let pr = match cfg.method {
                Method::SPP => cd[index]
                    .prefered_pseudorange()
                    .ok_or(Error::MissingPseudoRange)?,
                Method::CPP | Method::PPP => cd[index]
                    .pseudorange_combination()
                    .ok_or(Error::PseudoRangeCombination)?,
            };

            let (pr, frequency) = (pr.value, pr.carrier.frequency());

            // frequency dependent delay
            for delay in &cfg.int_delay {
                if delay.frequency == frequency {
                    models += delay.delay * SPEED_OF_LIGHT;
                }
            }

            /*
             * IONO + TROPO biases
             */
            let rtm = BiasRuntimeParams {
                t: cd[index].t,
                elevation,
                azimuth,
                frequency,
                apriori_geo,
            };

            /*
             * TROPO
             */
            if cfg.modeling.tropo_delay {
                if tropo_bias.needs_modeling() {
                    let bias = TroposphereBias::model(TropoModel::Niel, &rtm);
                    debug!("{} : modeled tropo delay {:.3E}[m]", cd[index].t, bias);
                    models += bias;
                    sv_input.tropo_bias = Bias::modeled(bias);
                } else if let Some(bias) = tropo_bias.bias(&rtm) {
                    debug!("{} : measured tropo delay {:.3E}[m]", cd[index].t, bias);
                    models += bias;
                    sv_input.tropo_bias = Bias::measured(bias);
                }
            }

            /*
             * IONO
             */
            if cfg.method == Method::SPP && cfg.modeling.iono_delay {
                if let Some(bias) = iono_bias.bias(&rtm) {
                    debug!(
                        "{} : modeled iono delay (f={:.3E}Hz) {:.3E}[m]",
                        cd[index].t, rtm.frequency, bias
                    );
                    models += bias;
                    sv_input.iono_bias = Bias::modeled(bias);
                }
            }

            y[i] = pr - rho - models;

            if i > 3 {
                g[(i, i)] = 1.0_f64;

                if cfg.method == Method::PPP && i > 3 {
                    let ph = cd[index]
                        .phase_combination()
                        .ok_or(Error::PseudoRangeCombination)?;

                    // TODO: conclude windup
                    let windup = 0.0_f64;
                    y[i] = ph.value - rho - models - windup;
                }
            }

            if i < cd.len() {
                sv.insert(cd[i].sv, sv_input);
            }
        }

        let w = cfg
            .solver
            .weight_matrix(sv.values().map(|sv| sv.elevation).collect());

        debug!("y: {} g: {}, w: {}", y, g, w);
        Ok(Self { y, g, w, sv })
    }
}

#[derive(Debug, Clone)]
pub(crate) struct Navigation {
    filter: Filter,
    pending: Output,
    filter_state: Option<FilterState>,
}

impl Navigation {
    pub fn new(filter: Filter) -> Self {
        Self {
            filter,
            filter_state: None,
            pending: Default::default(),
        }
    }
    pub fn resolve(&mut self, input: &Input) -> Result<Output, Error> {
        let out = self.filter.resolve(input, self.filter_state.clone())?;
        self.pending = out.clone();
        Ok(out)
    }
    pub fn validate(&mut self) {
        self.filter_state = Some(self.pending.state.clone());
    }
}
