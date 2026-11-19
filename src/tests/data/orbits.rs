use crate::{
    prelude::{Almanac, Duration, Epoch, Frame, Orbit, OrbitSource, SV},
    tests::data::{E01, E03, E05, E09, E13, E15, E24, E31},
};

use rstest::*;

use itertools::Itertools;

use std::{collections::HashMap, str::FromStr};

const MIN_ORBITS_PER_EPOCH: usize = 4 * 2;

#[fixture]
fn build_almanac() -> Almanac {
    use crate::tests::almanac;
    almanac()
}

#[fixture]
fn build_earth_frame() -> Frame {
    use crate::tests::earth_frame;
    earth_frame()
}

#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Eq, Hash)]
pub struct OrbitsDataKey {
    pub sv: SV,
    pub epoch: Epoch,
}

pub struct OrbitsData {
    pub map: HashMap<OrbitsDataKey, Orbit>,
}

impl OrbitsData {
    pub fn new(frame: Frame) -> Self {
        let t0_gpst = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();
        let t1_gpst = Epoch::from_str("2020-06-25T00:15:00 GPST").unwrap();
        let t2_gpst = Epoch::from_str("2020-06-25T00:30:00 GPST").unwrap();
        let t3_gpst = Epoch::from_str("2020-06-25T00:45:00 GPST").unwrap();
        let t4_gpst = Epoch::from_str("2020-06-25T01:00:00 GPST").unwrap();

        let map = HashMap::from_iter([
            (
                OrbitsDataKey {
                    sv: E01,
                    epoch: t0_gpst,
                },
                Orbit::from_position(-11562.163582, 14053.114306, 23345.128269, t0_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E03,
                    epoch: t0_gpst,
                },
                Orbit::from_position(4577.136069, -22995.974895, 18062.640686, t0_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E05,
                    epoch: t0_gpst,
                },
                Orbit::from_position(16577.017768, -4619.539763, 24092.494804, t0_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E09,
                    epoch: t0_gpst,
                },
                Orbit::from_position(18846.610510, 16144.830741, 16159.863309, t0_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E13,
                    epoch: t0_gpst,
                },
                Orbit::from_position(-15921.765341, -5400.108297, 24360.804625, t0_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E15,
                    epoch: t0_gpst,
                },
                Orbit::from_position(-409.147663, -21456.140629, 20391.202816, t0_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E24,
                    epoch: t0_gpst,
                },
                Orbit::from_position(26947.918853, 9320.740084, 7908.580547, t0_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E31,
                    epoch: t0_gpst,
                },
                Orbit::from_position(11195.440434, 16391.022663, 21968.198345, t0_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E01,
                    epoch: t1_gpst,
                },
                Orbit::from_position(-13618.625154, 13865.251337, 22325.739925, t1_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E03,
                    epoch: t1_gpst,
                },
                Orbit::from_position(5249.050928, -21427.275276, 19733.539343, t1_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E05,
                    epoch: t1_gpst,
                },
                Orbit::from_position(17451.377203, -2644.578813, 23770.958635, t1_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E09,
                    epoch: t1_gpst,
                },
                Orbit::from_position(19437.621892, 17372.955840, 14048.13756, t1_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E13,
                    epoch: t1_gpst,
                },
                Orbit::from_position(-14620.483100, -7198.531437, 24709.457156, t1_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E15,
                    epoch: t1_gpst,
                },
                Orbit::from_position(549.559454, -22943.716677, 18698.434257, t1_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E24,
                    epoch: t1_gpst,
                },
                Orbit::from_position(26132.488800, 9133.294721, 10452.358017, t1_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E31,
                    epoch: t1_gpst,
                },
                Orbit::from_position(9177.338648, 16128.279322, 23070.547459, t1_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E01,
                    epoch: t2_gpst,
                },
                Orbit::from_position(-15578.906571, 13828.422191, 21028.690065, t2_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E03,
                    epoch: t2_gpst,
                },
                Orbit::from_position(6089.334963, -19783.358422, 21158.979160, t2_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E05,
                    epoch: t2_gpst,
                },
                Orbit::from_position(18436.719493, -770.324223, 23154.020400, t2_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E09,
                    epoch: t2_gpst,
                },
                Orbit::from_position(20017.686424, 18382.862213, 11761.922747, t2_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E13,
                    epoch: t2_gpst,
                },
                Orbit::from_position(-13443.345950, -9101.204022, 24750.808535, t2_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E15,
                    epoch: t2_gpst,
                },
                Orbit::from_position(1313.090113, -24357.643799, 16773.179621, t2_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E24,
                    epoch: t2_gpst,
                },
                Orbit::from_position(25082.310095, 8993.926987, 12865.989009, t2_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E31,
                    epoch: t2_gpst,
                },
                Orbit::from_position(7059.263844, 16006.140471, 23886.192153, t2_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E01,
                    epoch: t3_gpst,
                },
                Orbit::from_position(-17408.137167, 13927.983030, 19470.096085, t3_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E03,
                    epoch: t3_gpst,
                },
                Orbit::from_position(7098.698927, -18099.833762, 22321.267348, t3_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E05,
                    epoch: t3_gpst,
                },
                Orbit::from_position(19509.930356, 974.904114, 22249.348845, t3_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E09,
                    epoch: t3_gpst,
                },
                Orbit::from_position(20553.093095, 19169.983128, 9329.608871, t3_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E13,
                    epoch: t3_gpst,
                },
                Orbit::from_position(-12411.471925, -11074.923521, 24484.357047, t3_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E15,
                    epoch: t3_gpst,
                },
                Orbit::from_position(1887.325160, -25661.135185, 14639.374285, t3_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E24,
                    epoch: t3_gpst,
                },
                Orbit::from_position(23814.294414, 8933.516592, 15119.417055, t3_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E31,
                    epoch: t3_gpst,
                },
                Orbit::from_position(4877.703417, 16036.340787, 24404.993554, t3_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E01,
                    epoch: t4_gpst,
                },
                Orbit::from_position(-19074.795786, 14143.814957, 17669.329791, t4_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E03,
                    epoch: t4_gpst,
                },
                Orbit::from_position(8271.797755, -16412.322261, 23205.987823, t4_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E05,
                    epoch: t4_gpst,
                },
                Orbit::from_position(20643.741422, 2567.147607, 21068.184349, t4_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E09,
                    epoch: t4_gpst,
                },
                Orbit::from_position(21010.298836, 19735.979210, 6781.401622, t4_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E13,
                    epoch: t4_gpst,
                },
                Orbit::from_position(-11540.233821, -13083.511580, 23913.427261, t4_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E15,
                    epoch: t4_gpst,
                },
                Orbit::from_position(2284.385592, -26819.485051, 12323.546523, t4_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E24,
                    epoch: t4_gpst,
                },
                Orbit::from_position(22350.983090, 8979.707681, 17184.581953, t4_gpst, frame),
            ),
            (
                OrbitsDataKey {
                    sv: E31,
                    epoch: t4_gpst,
                },
                Orbit::from_position(2670.799462, 16224.467677, 24620.497432, t4_gpst, frame),
            ),
        ]);

        Self { map }
    }
}

impl OrbitSource for OrbitsData {
    /// The solver requests the orbital state at transmission time (about 70 ms
    /// before the sampling epoch): the frozen SP3 samples (15') are interpolated
    /// (Lagrange, all samples of this satellite) at the requested instant.
    /// Returning the sample of the reception epoch would misplace each satellite
    /// by ~300 m along its track.
    fn state_at(&self, epoch: Epoch, sv: SV, frame: Frame) -> Option<Orbit> {
        let samples = self
            .map
            .iter()
            .filter(|(k, _)| k.sv == sv)
            .map(|(k, v)| {
                let pos_vel_km = v.to_cartesian_pos_vel();
                (
                    (k.epoch - epoch).to_seconds(),
                    [pos_vel_km[0], pos_vel_km[1], pos_vel_km[2]],
                )
            })
            .collect::<Vec<_>>();

        let (x_km, y_km, z_km) = lagrange_interpolation(&samples)?;

        Some(Orbit::from_position(x_km, y_km, z_km, epoch, frame))
    }
}

/// Maximal extrapolation (s): first epoch transmission times precede the first sample.
pub(crate) const MAX_EXTRAPOLATION_S: f64 = 300.0;

/// Lagrange interpolation at t = 0, of (t_i, [x_i, y_i, z_i]) samples.
/// Returns None if less than 2 samples, or if t = 0 is more than
/// [MAX_EXTRAPOLATION_S] outside of the sampling span.
pub(crate) fn lagrange_interpolation(samples: &[(f64, [f64; 3])]) -> Option<(f64, f64, f64)> {
    if samples.len() < 2 {
        return None;
    }

    let t_min = samples
        .iter()
        .map(|(t, _)| *t)
        .fold(f64::INFINITY, f64::min);
    let t_max = samples
        .iter()
        .map(|(t, _)| *t)
        .fold(f64::NEG_INFINITY, f64::max);

    if t_min > MAX_EXTRAPOLATION_S || t_max < -MAX_EXTRAPOLATION_S {
        return None;
    }

    let mut out = [0.0; 3];

    for (i, (t_i, p_i)) in samples.iter().enumerate() {
        let mut l_i = 1.0;

        for (j, (t_j, _)) in samples.iter().enumerate() {
            if i != j {
                l_i *= (0.0 - t_j) / (t_i - t_j);
            }
        }

        for k in 0..3 {
            out[k] += l_i * p_i[k];
        }
    }

    Some((out[0], out[1], out[2]))
}

#[test]
fn frozen_orbits_interpolation() {
    let earth_frame = build_earth_frame();
    let data = OrbitsData::new(earth_frame);

    let t1 = Epoch::from_str("2020-06-25T00:15:00 GPST").unwrap();

    // on a sample: exact
    let sample = data.map[&OrbitsDataKey { sv: E01, epoch: t1 }].to_cartesian_pos_vel();
    let interp = data
        .state_at(t1, E01, earth_frame)
        .unwrap()
        .to_cartesian_pos_vel();

    for i in 0..3 {
        assert!((sample[i] - interp[i]).abs() < 1E-6);
    }

    // 75 ms earlier: the satellite moved by v.dt (earth fixed velocity: ~200 to 300 m)
    let tx = t1 - Duration::from_milliseconds(75.0);
    let moved = data
        .state_at(tx, E01, earth_frame)
        .unwrap()
        .to_cartesian_pos_vel();

    let d_m = ((moved[0] - sample[0]).powi(2)
        + (moved[1] - sample[1]).powi(2)
        + (moved[2] - sample[2]).powi(2))
    .sqrt()
        * 1.0E3;
    assert!(d_m > 150.0 && d_m < 350.0, "satellite moved by {d_m} m");

    // no extrapolation
    let outside = Epoch::from_str("2020-06-25T02:00:00 GPST").unwrap();
    assert!(data.state_at(outside, E01, earth_frame).is_none());
}

#[test]
fn verify_min_test_orbits_per_epoch() {
    let earth_frame = build_earth_frame();

    let data = OrbitsData::new(earth_frame);

    for epoch in data.map.keys().map(|k| k.epoch).unique() {
        let collected = data
            .map
            .iter()
            .filter(|(k, _)| k.epoch == epoch)
            .collect::<Vec<_>>();

        assert!(
            collected.len() >= MIN_ORBITS_PER_EPOCH,
            "not enough valid test orbits @ {epoch}"
        );
    }
}
