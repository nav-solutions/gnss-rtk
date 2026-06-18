use log::error;

use crate::{
    candidate::differences::Differences,
    navigation::{sv::SVContribution, vector::VectorContribution},
    prelude::{Candidate, Config, Epoch, Error, Method, Vector3},
};

impl Candidate {
    /// Measurement vector contribution.
    pub(crate) fn rtk_vector_contribution(
        &self,
        epoch: Epoch,
        two_rows: bool,
        cfg: &Config,
        double_diffs: &Differences,
        contribution: &mut SVContribution,
    ) -> Result<VectorContribution, Error> {
        let dd = double_diffs
            .difference(self.sv)
            .ok_or(Error::RtkDDPostfitMissing)?;

        let mut vec = VectorContribution::default();

        // row #1
        match cfg.method {
            Method::SPP => {
                if let Some((_, code)) = dd.code {
                    vec.row_1 = code;
                } else {
                    error!("{}({}) - missing pseudo range", epoch, self.sv);
                    return Err(Error::MissingPseudoRange);
                }
            },
            _ => {
                if let Some((_, code)) = dd.code_if {
                    vec.row_1 = code;
                } else {
                    error!("{}({}) - missing pseudo range", epoch, self.sv);
                    return Err(Error::MissingPseudoRange);
                }
            },
        }

        // row #1
        if !two_rows && cfg.method == Method::PPP {
            if let Some(phase_if) = dd.phase_if(self.sv) {
                vec.row_1 = phase_if;
            } else {
                error!("{}({}) - missing phase data", epoch, self.sv);
                return Err(Error::MissingPhaseRange);
            }
        }

        // row #2 (special case)
        if two_rows && cfg.method == Method::PPP {
            // special case
            if let Some((_, _, phase_if)) = dd.phase_if {
                vec.row_2 = phase_if;
            } else {
                error!("{}({}) - missing phase data", epoch, self.sv);
                return Err(Error::MissingPhaseRange);
            }
        }

        // Measurement weighting, keyed off this satellite's elevation.
        // NOTE: first-order model applied per double-difference row; it does
        // not yet account for DD cross-correlation through the shared pivot.
        let elevation_deg = match self.attitude() {
            Some((elev_deg, azim_deg)) => {
                contribution.elevation_deg = elev_deg;
                contribution.azimuth_deg = azim_deg;
                Some(elev_deg)
            },
            None => None,
        };

        if !two_rows && cfg.method == Method::PPP {
            // row #1 carries phase in this case
            vec.sigma_1 = cfg.solver.phase_sigma_m(elevation_deg);
        } else {
            vec.sigma_1 = cfg.solver.code_sigma_m(elevation_deg);
        }

        if two_rows && cfg.method == Method::PPP {
            vec.sigma_2 = cfg.solver.phase_sigma_m(elevation_deg);
        }

        Ok(vec)
    }

    /// Matrix contribution.
    ///
    /// ## Input
    ///  - i: matrix row
    ///  - cfg: [Config] preset
    ///  - x0_y0_z0: position coordinates as ECEF (m)
    pub(crate) fn rtk_matrix_contribution(
        &self,
        x0_y0_z0_m: Vector3<f64>,
        pivot_position_ecef_m: (f64, f64, f64),
    ) -> (f64, f64, f64) {
        let (x0_m, y0_m, z0_m) = (x0_y0_z0_m[0], x0_y0_z0_m[1], x0_y0_z0_m[2]);
        let (pivot_x_m, pivot_y_m, pivot_z_m) = pivot_position_ecef_m;

        let orbit = self.orbit.unwrap_or_else(|| {
            panic!(
                "internal error: {}({}) state not fully resolved!",
                self.epoch, self.sv
            );
        });

        let pos_vel_m = orbit.to_cartesian_pos_vel() * 1.0E3;
        let (sv_x_m, sv_y_m, sv_z_m) = (pos_vel_m[0], pos_vel_m[1], pos_vel_m[2]);

        let rho_i =
            ((sv_x_m - x0_m).powi(2) + (sv_y_m - y0_m).powi(2) + (sv_z_m - z0_m).powi(2)).sqrt();

        let e_i = (
            (sv_x_m - x0_m) / rho_i,
            (sv_y_m - y0_m) / rho_i,
            (sv_z_m - z0_m) / rho_i,
        );

        let rho_j =
            ((pivot_x_m - x0_m).powi(2) + (pivot_y_m - y0_m).powi(2) + (pivot_z_m - z0_m).powi(2))
                .sqrt();

        let e_j = (
            (pivot_x_m - x0_m) / rho_j,
            (pivot_y_m - y0_m) / rho_j,
            (pivot_z_m - z0_m) / rho_j,
        );

        (e_j.0 - e_i.0, e_j.1 - e_i.1, e_j.2 - e_i.2)
    }
}

#[cfg(test)]
mod test {
    use crate::{
        prelude::{Epoch, Frame, Orbit},
        tests::{CandidatesBuilder, E05, ROVER_REFERENCE_COORDS_ECEF_M},
    };

    use nalgebra::Vector3;
    use rstest::*;

    use std::str::FromStr;

    #[fixture]
    fn build_earth_frame() -> Frame {
        use crate::tests::earth_frame;
        earth_frame()
    }

    #[test]
    fn rtk_matrix_contribution() {
        let earth_frame = build_earth_frame();
        let t0 = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

        let pivot_position_ecef_km = (16577.017768, -4619.539763, 24092.494804);

        let pivot_position_ecef_m = (
            pivot_position_ecef_km.0 * 1000.0,
            pivot_position_ecef_km.1 * 1000.0,
            pivot_position_ecef_km.2 * 1000.0,
        );

        let e01_position_ecef_km = (-11562.163582, 14053.114306, 23345.128269);

        let e01_position_ecef_m = (
            e01_position_ecef_km.0 * 1000.0,
            e01_position_ecef_km.1 * 1000.0,
            e01_position_ecef_km.2 * 1000.0,
        );

        let mut rover = CandidatesBuilder::build_rover_sv_at(E05, t0);

        rover.orbit = Some(Orbit::from_position(
            e01_position_ecef_km.0,
            e01_position_ecef_km.1,
            e01_position_ecef_km.2,
            t0,
            earth_frame,
        ));

        let x0_y0_z0_m = Vector3::new(
            ROVER_REFERENCE_COORDS_ECEF_M.0,
            ROVER_REFERENCE_COORDS_ECEF_M.1,
            ROVER_REFERENCE_COORDS_ECEF_M.2,
        );

        let (dx, dy, dz) = rover.rtk_matrix_contribution(x0_y0_z0_m, pivot_position_ecef_m);

        let rho_i = ((e01_position_ecef_m.0 - ROVER_REFERENCE_COORDS_ECEF_M.0).powi(2)
            + (e01_position_ecef_m.1 - ROVER_REFERENCE_COORDS_ECEF_M.1).powi(2)
            + (e01_position_ecef_m.2 - ROVER_REFERENCE_COORDS_ECEF_M.2).powi(2))
        .sqrt();

        let rho_j = ((pivot_position_ecef_m.0 - ROVER_REFERENCE_COORDS_ECEF_M.0).powi(2)
            + (pivot_position_ecef_m.1 - ROVER_REFERENCE_COORDS_ECEF_M.1).powi(2)
            + (pivot_position_ecef_m.2 - ROVER_REFERENCE_COORDS_ECEF_M.2).powi(2))
        .sqrt();

        let e_i = (
            (e01_position_ecef_m.0 - ROVER_REFERENCE_COORDS_ECEF_M.0) / rho_i,
            (e01_position_ecef_m.1 - ROVER_REFERENCE_COORDS_ECEF_M.1) / rho_i,
            (e01_position_ecef_m.2 - ROVER_REFERENCE_COORDS_ECEF_M.2) / rho_i,
        );

        let e_j = (
            (pivot_position_ecef_m.0 - ROVER_REFERENCE_COORDS_ECEF_M.0) / rho_j,
            (pivot_position_ecef_m.1 - ROVER_REFERENCE_COORDS_ECEF_M.1) / rho_j,
            (pivot_position_ecef_m.2 - ROVER_REFERENCE_COORDS_ECEF_M.2) / rho_j,
        );

        assert!((dx - (e_j.0 - e_i.0)).abs() < 1E-6,);
        assert!((dy - (e_j.1 - e_i.1)).abs() < 1E-6,);
        assert!((dz - (e_j.2 - e_i.2)).abs() < 1E-6,);
    }
}
