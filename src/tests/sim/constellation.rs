//! Synthetic constellations (Walker circular orbits), expressed in ECEF.
//!
//! The orbits are not meant to be realistic at the meter level (no J2, no
//! eccentricity), they only need to be smooth, physically plausible, and known
//! exactly by the simulator. The solver only sees them through [OrbitSource].
use nalgebra::Vector3;

use crate::{
    constants::EARTH_ANGULAR_VEL_RAD,
    prelude::{Constellation, SV},
};

/// Earth gravitational parameter (m³.s⁻²)
const GM_M3_S2: f64 = 3.986_004_418E14;

/// One circular orbit
#[derive(Debug, Clone, Copy)]
pub struct CircularOrbit {
    pub sv: SV,
    /// Semi major axis (m)
    pub a_m: f64,
    /// Inclination (rad)
    pub inc_rad: f64,
    /// Right ascension of ascending node at t=0 (rad)
    pub raan_rad: f64,
    /// Argument of latitude at t=0 (rad)
    pub u0_rad: f64,
}

impl CircularOrbit {
    /// Mean motion (rad/s)
    pub fn mean_motion(&self) -> f64 {
        (GM_M3_S2 / self.a_m.powi(3)).sqrt()
    }

    /// Inertial position (m), `t_s` seconds after scenario start.
    pub fn eci_position_m(&self, t_s: f64) -> Vector3<f64> {
        let u = self.u0_rad + self.mean_motion() * t_s;
        let (su, cu) = u.sin_cos();
        let (so, co) = self.raan_rad.sin_cos();
        let (si, ci) = self.inc_rad.sin_cos();

        Vector3::new(
            self.a_m * (cu * co - su * ci * so),
            self.a_m * (cu * so + su * ci * co),
            self.a_m * (su * si),
        )
    }

    /// Earth fixed position (m), `t_s` seconds after scenario start.
    pub fn ecef_position_m(&self, t_s: f64) -> Vector3<f64> {
        let r = self.eci_position_m(t_s);
        let theta = EARTH_ANGULAR_VEL_RAD * t_s;
        let (s, c) = theta.sin_cos();
        Vector3::new(c * r[0] + s * r[1], -s * r[0] + c * r[1], r[2])
    }
}

/// Supported synthetic constellations
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum SimConstellation {
    /// Galileo like Walker 24/3/1, a=29600 km, i=56°
    Galileo,
    /// GPS like 24 satellites on 6 planes, a=26560 km, i=55°
    GPS,
}

impl SimConstellation {
    pub fn orbits(&self) -> Vec<CircularOrbit> {
        match self {
            Self::Galileo => walker(Constellation::Galileo, 29_599_801.0, 56.0, 24, 3, 1, 0.3),
            Self::GPS => walker(Constellation::GPS, 26_559_710.0, 55.0, 24, 6, 1, 1.1),
        }
    }
}

/// Walker delta pattern t/p/f
fn walker(
    constellation: Constellation,
    a_m: f64,
    inc_deg: f64,
    t: usize,
    p: usize,
    f: usize,
    raan_offset_rad: f64,
) -> Vec<CircularOrbit> {
    let s = t / p;
    let two_pi = 2.0 * std::f64::consts::PI;

    let mut out = Vec::with_capacity(t);

    for plane in 0..p {
        for k in 0..s {
            let prn = (plane * s + k + 1) as u8;
            let raan = raan_offset_rad + two_pi * plane as f64 / p as f64;
            let u0 = two_pi * k as f64 / s as f64 + two_pi * (f * plane) as f64 / t as f64;

            out.push(CircularOrbit {
                sv: SV::new(constellation, prn),
                a_m,
                inc_rad: inc_deg.to_radians(),
                raan_rad: raan,
                u0_rad: u0,
            });
        }
    }

    out
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn orbital_radius_and_period() {
        for c in [SimConstellation::Galileo, SimConstellation::GPS] {
            let orbits = c.orbits();
            assert_eq!(orbits.len(), 24);
            for orb in orbits {
                for t in [0.0, 1000.0, 43200.0] {
                    let r = orb.ecef_position_m(t).norm();
                    assert!((r - orb.a_m).abs() < 1E-3);
                }
                let period_h = 2.0 * std::f64::consts::PI / orb.mean_motion() / 3600.0;
                assert!(period_h > 11.0 && period_h < 15.0, "period={period_h}h");
            }
        }
    }

    #[test]
    fn satellite_speed_is_plausible() {
        let orb = SimConstellation::GPS.orbits()[0];
        let v = (orb.eci_position_m(1.0) - orb.eci_position_m(0.0)).norm();
        assert!((v - 3874.0).abs() < 20.0, "v={v} m/s");
    }
}
