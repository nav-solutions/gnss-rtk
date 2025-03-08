use nalgebra::{base::dimension::U4, ArrayStorage, Matrix, Matrix3, Matrix4};

use crate::navigation::state::State;

/// [Navigation] filter [DilutionOfPrecision]
#[derive(Clone, Default, Copy)]
pub(crate) struct DilutionOfPrecision {
    /// Geometric DOP
    pub gdop: f64,
    /// Horizontal DOP
    pub hdop: f64,
    /// Vertical DOP
    pub vdop: f64,
    /// Temporal DOP
    pub tdop: f64,
}

impl DilutionOfPrecision {
    pub(crate) fn q_enu(h: Matrix4<f64>, lat_rad: f64, lon_rad: f64) -> Matrix3<f64> {
        let r = Matrix3::<f64>::new(
            -lon_rad.sin(),
            -lon_rad.cos() * lat_rad.sin(),
            lat_rad.cos() * lon_rad.cos(),
            lon_rad.cos(),
            -lat_rad.sin() * lon_rad.sin(),
            lat_rad.cos() * lon_rad.sin(),
            0.0_f64,
            lat_rad.cos(),
            lon_rad.sin(),
        );

        let q_3 = Matrix3::<f64>::new(
            h[(0, 0)],
            h[(0, 1)],
            h[(0, 2)],
            h[(1, 0)],
            h[(1, 1)],
            h[(1, 2)],
            h[(2, 0)],
            h[(2, 1)],
            h[(2, 2)],
        );

        r.clone().transpose() * q_3 * r
    }

    /// Creates new [DillutionOfPrecision] from matrix
    pub fn new(state: &State, g: Matrix<f64, U4, U4, ArrayStorage<f64, 4, 4>>) -> Self {
        let (lat_rad, long_rad) = (
            state.lat_long_alt_deg_deg_km.0.to_radians(),
            state.lat_long_alt_deg_deg_km.1.to_radians(),
        );

        let q_enu = Self::q_enu(g, lat_rad, long_rad);

        Self {
            gdop: g.trace().sqrt(),
            tdop: g[(3, 3)].sqrt(),
            vdop: q_enu[(2, 2)].sqrt(),
            hdop: (q_enu[(0, 0)] + q_enu[(1, 1)]).sqrt(),
        }
    }
}
