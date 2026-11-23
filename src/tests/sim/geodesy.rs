//! Minimal WGS84 geodesy helpers used by the simulator (independent of the solver).
use nalgebra::{Matrix3, Vector3};

pub const WGS84_A_M: f64 = 6_378_137.0;
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;

fn e2() -> f64 {
    WGS84_F * (2.0 - WGS84_F)
}

/// Geodetic (lat rad, lon rad, height m) to ECEF (m).
pub fn geodetic_to_ecef(lat_rad: f64, lon_rad: f64, h_m: f64) -> Vector3<f64> {
    let (slat, clat) = lat_rad.sin_cos();
    let (slon, clon) = lon_rad.sin_cos();
    let n = WGS84_A_M / (1.0 - e2() * slat * slat).sqrt();

    Vector3::new(
        (n + h_m) * clat * clon,
        (n + h_m) * clat * slon,
        (n * (1.0 - e2()) + h_m) * slat,
    )
}

/// ECEF (m) to geodetic (lat rad, lon rad, height m), iterative.
pub fn ecef_to_geodetic(r: &Vector3<f64>) -> (f64, f64, f64) {
    let (x, y, z) = (r[0], r[1], r[2]);
    let lon = y.atan2(x);
    let p = (x * x + y * y).sqrt();

    let mut lat = z.atan2(p * (1.0 - e2()));
    let mut h = 0.0;

    for _ in 0..10 {
        let slat = lat.sin();
        let n = WGS84_A_M / (1.0 - e2() * slat * slat).sqrt();
        h = p / lat.cos() - n;
        lat = z.atan2(p * (1.0 - e2() * n / (n + h)));
    }

    (lat, lon, h)
}

/// ECEF -> ENU rotation matrix at given geodetic location.
pub fn ecef_to_enu_matrix(lat_rad: f64, lon_rad: f64) -> Matrix3<f64> {
    let (slat, clat) = lat_rad.sin_cos();
    let (slon, clon) = lon_rad.sin_cos();

    Matrix3::new(
        -slon,
        clon,
        0.0,
        -slat * clon,
        -slat * slon,
        clat,
        clat * clon,
        clat * slon,
        slat,
    )
}

/// Expresses an ECEF difference vector in the local ENU frame of `reference`.
pub fn ecef_delta_to_enu(reference: &Vector3<f64>, delta: &Vector3<f64>) -> Vector3<f64> {
    let (lat, lon, _) = ecef_to_geodetic(reference);
    ecef_to_enu_matrix(lat, lon) * delta
}

/// (elevation, azimuth) in degrees of `target` seen from `observer` (both ECEF m).
pub fn elevation_azimuth_deg(observer: &Vector3<f64>, target: &Vector3<f64>) -> (f64, f64) {
    let enu = ecef_delta_to_enu(observer, &(target - observer));
    let horiz = (enu[0] * enu[0] + enu[1] * enu[1]).sqrt();
    let el = enu[2].atan2(horiz).to_degrees();
    let mut az = enu[0].atan2(enu[1]).to_degrees();
    if az < 0.0 {
        az += 360.0;
    }
    (el, az)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn geodetic_roundtrip() {
        for (lat, lon, h) in [
            (55.5_f64, 8.45_f64, 50.0),
            (-33.9, 151.2, 10.0),
            (0.0, 0.0, 0.0),
            (89.0, -120.0, 3000.0),
        ] {
            let r = geodetic_to_ecef(lat.to_radians(), lon.to_radians(), h);
            let (la, lo, hh) = ecef_to_geodetic(&r);
            assert!((la.to_degrees() - lat).abs() < 1E-9);
            assert!((lo.to_degrees() - lon).abs() < 1E-9);
            assert!((hh - h).abs() < 1E-4);
        }
    }

    #[test]
    fn zenith_elevation() {
        let rx = geodetic_to_ecef(0.7, 0.2, 100.0);
        let up = rx + rx.normalize() * 20_000_000.0;
        let (el, _) = elevation_azimuth_deg(&rx, &up);
        assert!(el > 89.5);
    }
}
