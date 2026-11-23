use crate::prelude::{BiasRuntime, Duration};

#[cfg(doc)]
use crate::prelude::Method;

/// Satellite clock correction.
///
/// [Duration] has a 1 ns resolution (0.3 m): precise products (SP3 / CLK, 1E-12 s
/// resolution) must be provided with the `*_s` constructors, which keep the sub
/// nanosecond part. A 1 ns quantization is harmless for code navigation, but its
/// ±0.15 m random error per satellite and per epoch ruins phase navigation (PPP).
#[derive(Default, Debug, Copy, Clone)]
pub struct SatelliteClockCorrection {
    /// Correction to associated timescale, expressed as [Duration]
    /// (truncated to the nanosecond, see [Self::seconds] for the exact value).
    pub duration: Duration,

    /// True when relativistic effect has not been corrected.
    pub(crate) needs_relativistic_correction: bool,

    /// Sub nanosecond remainder (s): exact value = duration + remainder
    remainder_s: f64,
}

impl SatelliteClockCorrection {
    /// Define a new [SatelliteClockCorrection] that already integrates relativistic corrections
    pub fn with_relativistic_correction(duration: Duration) -> Self {
        Self {
            duration,
            needs_relativistic_correction: false,
            remainder_s: 0.0,
        }
    }

    /// Define a new [SatelliteClockCorrection] that does not integrate relativistic corrections
    pub fn without_relativistic_correction(duration: Duration) -> Self {
        Self {
            duration,
            needs_relativistic_correction: true,
            remainder_s: 0.0,
        }
    }

    /// Define a new [SatelliteClockCorrection] that already integrates relativistic
    /// corrections, from a value in seconds, without loss of precision.
    pub fn with_relativistic_correction_s(seconds: f64) -> Self {
        let mut s = Self::with_relativistic_correction(Duration::from_seconds(seconds));
        s.remainder_s = seconds - s.duration.to_seconds();
        s
    }

    /// Define a new [SatelliteClockCorrection] that does not integrate relativistic
    /// corrections, from a value in seconds, without loss of precision.
    pub fn without_relativistic_correction_s(seconds: f64) -> Self {
        let mut s = Self::without_relativistic_correction(Duration::from_seconds(seconds));
        s.remainder_s = seconds - s.duration.to_seconds();
        s
    }

    /// Exact clock correction, in seconds (sub nanosecond part included).
    pub fn seconds(&self) -> f64 {
        self.duration.to_seconds() + self.remainder_s
    }
}

/// [SpacebornBias] must be implemented to determine all on-board biases.
/// Basic navigation only requires [SpacebornBias::satellite_clock_bias],
/// but precise navigation requires all of it.
pub trait SpacebornBias {
    /// Provide the [SatelliteClockCorrection] for requested satellite,
    /// as epoch of observation.
    /// If this is not known, simply return [Default::default()], but
    /// it will dramatically impact the accuracy of any aboslute solution (>100km).
    /// This value is disregarded when sv_clock_bias is not being modeled.
    /// This value does not impact the accuracy of RTK solutions.
    fn clock_bias(&self, rtm: &BiasRuntime) -> SatelliteClockCorrection;

    /// Provide the Satellite specific on-board group delay,
    /// referece to the L1 frequency.
    /// This value will be disregarded when sv_total_group_delay
    /// is not being modeled. If this value is not known,
    /// simply return Duration::ZERO.
    /// This is mandatory when sv_total_group_delay is being modeled,
    /// otherwise this [SV] will get dropped.
    fn group_delay(&self, rtm: &BiasRuntime) -> Duration;

    /// Provide the Satellite MW internal bias for requested [SV].
    /// This is mandatory in [Method::PPP] technique is being used,
    /// otherwise this [SV] will get dropped.
    fn mw_bias(&self, rtm: &BiasRuntime) -> f64;
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn sub_nanosecond_clock_corrections() {
        // GRG clock product value (1E-12 s resolution)
        let value_s = -1.533_345_678_912E-5;

        // Duration: truncated to the nanosecond (up to 0.3 m)
        let legacy =
            SatelliteClockCorrection::with_relativistic_correction(Duration::from_seconds(value_s));
        assert!((legacy.seconds() - value_s).abs() < 1.0E-9);
        assert!((legacy.seconds() - value_s).abs() > 1.0E-12);

        // exact
        for clock in [
            SatelliteClockCorrection::with_relativistic_correction_s(value_s),
            SatelliteClockCorrection::without_relativistic_correction_s(value_s),
        ] {
            assert!((clock.seconds() - value_s).abs() < 1.0E-15, "{clock:?}");
            assert!((clock.duration.to_seconds() - value_s).abs() < 1.0E-9);
        }

        assert!(
            !SatelliteClockCorrection::with_relativistic_correction_s(value_s)
                .needs_relativistic_correction
        );
        assert!(
            SatelliteClockCorrection::without_relativistic_correction_s(value_s)
                .needs_relativistic_correction
        );

        assert_eq!(SatelliteClockCorrection::default().seconds(), 0.0);
    }
}
