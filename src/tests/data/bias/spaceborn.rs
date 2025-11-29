use crate::{
    prelude::{BiasRuntime, Duration, Epoch, SV, SatelliteClockCorrection, SpacebornBias},
    tests::{E01, E03, E05, E09, E13, E15, E24, E31},
};

use std::collections::HashMap;
use std::str::FromStr;

#[derive(Clone, Copy)]
struct SVBias {
    pub tgd: Duration,
    pub mw_bias: f64,
    pub clock_bias: SatelliteClockCorrection,
}

pub struct TestSpacebornBiases {
    sv_biases: HashMap<(SV, Epoch), SVBias>,
}

impl TestSpacebornBiases {
    pub fn build() -> Self {
        let dataset = [
            (
                (E01, "2020-06-25T00:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_seconds(-0.884707516318E-03),
                    ),
                },
            ),
            (
                (E03, "2020-06-25T00:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_seconds(-0.313499770596E-03),
                    ),
                },
            ),
            (
                (E05, "2020-06-25T00:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_seconds(-0.368776159133E-03),
                    ),
                },
            ),
            (
                (E09, "2020-06-25T00:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_seconds(0.601769391412E-02),
                    ),
                },
            ),
            (
                (E13, "2020-06-25T00:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_seconds(0.401847043889E-03),
                    ),
                },
            ),
            (
                (E15, "2020-06-25T00:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_seconds(0.862331400195E-03),
                    ),
                },
            ),
            (
                (E24, "2020-06-25T00:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_seconds(0.538503520147E-02),
                    ),
                },
            ),
            (
                (E31, "2020-06-25T00:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_seconds(-0.472987972013E-03),
                    ),
                },
            ),
            (
                (E01, "2020-06-25T00:15:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-884.714669),
                    ),
                },
            ),
            (
                (E03, "2020-06-25T00:15:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-313.503498),
                    ),
                },
            ),
            (
                (E05, "2020-06-25T00:15:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-368.773276),
                    ),
                },
            ),
            (
                (E09, "2020-06-25T00:15:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(6017.682886),
                    ),
                },
            ),
            (
                (E13, "2020-06-25T00:15:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(401.847317),
                    ),
                },
            ),
            (
                (E15, "2020-06-25T00:15:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(862.330144),
                    ),
                },
            ),
            (
                (E24, "2020-06-25T00:15:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(5385.017319),
                    ),
                },
            ),
            (
                (E31, "2020-06-25T00:15:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-472.988096),
                    ),
                },
            ),
            (
                (E01, "2020-06-25T00:30:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-884.721810),
                    ),
                },
            ),
            (
                (E03, "2020-06-25T00:30:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-313.507266),
                    ),
                },
            ),
            (
                (E05, "2020-06-25T00:30:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-368.770391),
                    ),
                },
            ),
            (
                (E09, "2020-06-25T00:30:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(6017.671834),
                    ),
                },
            ),
            (
                (E13, "2020-06-25T00:30:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(401.847606),
                    ),
                },
            ),
            (
                (E15, "2020-06-25T00:30:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(862.328904),
                    ),
                },
            ),
            (
                (E24, "2020-06-25T00:30:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(5384.999417),
                    ),
                },
            ),
            (
                (E31, "2020-06-25T00:30:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-472.988276),
                    ),
                },
            ),
            (
                (E01, "2020-06-25T00:45:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-884.728967),
                    ),
                },
            ),
            (
                (E03, "2020-06-25T00:45:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-313.511002),
                    ),
                },
            ),
            (
                (E05, "2020-06-25T00:45:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-368.767474),
                    ),
                },
            ),
            (
                (E09, "2020-06-25T00:45:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(6017.660829),
                    ),
                },
            ),
            (
                (E13, "2020-06-25T00:45:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(401.847841),
                    ),
                },
            ),
            (
                (E15, "2020-06-25T00:45:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(862.327697),
                    ),
                },
            ),
            (
                (E24, "2020-06-25T00:45:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(5384.981517),
                    ),
                },
            ),
            (
                (E31, "2020-06-25T00:45:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-472.988448),
                    ),
                },
            ),
            (
                (E03, "2020-06-25T01:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-313.514700),
                    ),
                },
            ),
            (
                (E05, "2020-06-25T01:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-368.764558),
                    ),
                },
            ),
            (
                (E09, "2020-06-25T01:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(6017.649787),
                    ),
                },
            ),
            (
                (E13, "2020-06-25T01:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(401.848117),
                    ),
                },
            ),
            (
                (E15, "2020-06-25T01:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(862.326442),
                    ),
                },
            ),
            (
                (E24, "2020-06-25T01:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(5384.963584),
                    ),
                },
            ),
            (
                (E31, "2020-06-25T01:00:00 GPST"),
                SVBias {
                    mw_bias: 0.0,
                    tgd: Duration::ZERO,
                    clock_bias: SatelliteClockCorrection::without_relativistic_correction(
                        Duration::from_microseconds(-472.988574),
                    ),
                },
            ),
        ];

        Self {
            sv_biases: HashMap::from_iter(dataset.into_iter().map(|((sv, t_str), bias)| {
                let t = Epoch::from_str(t_str).unwrap_or_else(|e| {
                    panic!("invalid epoch specs: {t_str} - {e}");
                });
                ((sv, t), bias)
            })),
        }
    }
}

impl SpacebornBias for TestSpacebornBiases {
    fn clock_bias(&self, rtm: &BiasRuntime) -> SatelliteClockCorrection {
        let sv_bias = self
            .sv_biases
            .iter()
            .filter_map(|(k, v)| {
                if k.0 == rtm.sv && k.1 == rtm.epoch {
                    Some(v)
                } else {
                    None
                }
            })
            .reduce(|k, _| k);

        if let Some(bias) = sv_bias {
            bias.clock_bias
        } else {
            Default::default()
        }
    }

    fn group_delay(&self, rtm: &BiasRuntime) -> Duration {
        let sv_bias = self
            .sv_biases
            .iter()
            .filter_map(|(k, v)| {
                if k.0 == rtm.sv && k.1 == rtm.epoch {
                    Some(v)
                } else {
                    None
                }
            })
            .reduce(|k, _| k);

        if let Some(bias) = sv_bias {
            bias.tgd
        } else {
            Duration::ZERO
        }
    }

    fn mw_bias(&self, rtm: &BiasRuntime) -> f64 {
        let sv_bias = self
            .sv_biases
            .iter()
            .filter_map(|(k, v)| {
                if k.0 == rtm.sv && k.1 == rtm.epoch {
                    Some(v)
                } else {
                    None
                }
            })
            .reduce(|k, _| k);

        if let Some(bias) = sv_bias {
            bias.mw_bias
        } else {
            0.0_f64
        }
    }
}

#[cfg(test)]
mod test {
    #[test]
    fn verify_space_bias() {}
}
