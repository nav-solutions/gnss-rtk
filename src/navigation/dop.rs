use nalgebra::{DMatrix, DimName, Matrix3, U3};

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
    pub(crate) fn q_enu(mat: &DMatrix<f64>, lat_rad: f64, lon_rad: f64) -> Matrix3<f64> {
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
            mat[(0, 0)],
            mat[(0, 1)],
            mat[(0, 2)],
            mat[(1, 0)],
            mat[(1, 1)],
            mat[(1, 2)],
            mat[(2, 0)],
            mat[(2, 1)],
            mat[(2, 2)],
        );

        r.clone().transpose() * q_3 * r
    }

    /// Creates new [DilutionOfPrecision].
    ///
    /// ## Inut
    /// - new [State]
    /// - g_g_t = (G * GT)⁻¹ matrix
    pub fn new(state: &State, g_gt_inv: DMatrix<f64>) -> Self {
        let (nrows, ncols) = (g_gt_inv.nrows(), g_gt_inv.ncols());

        assert!(nrows >= U3::DIM, "incorrect (G.G)⁻¹ dimensions");
        assert_eq!(nrows, ncols, "(G.G)⁻¹ is not square");

        let (lat_rad, long_rad) = (
            state.lat_long_alt_deg_deg_km.0.to_radians(),
            state.lat_long_alt_deg_deg_km.1.to_radians(),
        );

        let q_enu = Self::q_enu(&g_gt_inv, lat_rad, long_rad);

        Self {
            gdop: g_gt_inv.trace().sqrt(),
            tdop: if nrows > U3::DIM {
                g_gt_inv[(3, 3)].sqrt()
            } else {
                0.0
            },
            vdop: q_enu[(2, 2)].sqrt(),
            hdop: (q_enu[(0, 0)] + q_enu[(1, 1)]).sqrt(),
        }
    }
}
