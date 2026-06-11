//! Position solving candidate
use log::debug;

use crate::{
    constants::SPEED_OF_LIGHT_M_S,
    navigation::{sv::SVContribution, vector::VectorContribution},
    prelude::{Candidate, Config, Duration, Error, Method, Signal, Vector3},
};

impl Candidate {
    /// Measurement vector contribution.
    /// This will pass if
    /// - State has been previously resolved
    /// - Range estimate is available
    /// - Preset modeling are matched
    ///
    /// ## Input
    /// - t: [Epoch] of computation
    /// - cfg: [Config] preset
    /// - x0_y0_z0: current state (metric)
    /// - rx_lat_long_alt_ddeg_km: state as geodetic lat, long both
    /// in decimal degrees, and altitude above mean sea level (km)
    /// - contribution: mutable [SVContribution]
    /// ## Returns
    /// - b_i contribution, r_i contribution, dr: relativistic path range
    pub(crate) fn ppp_vector_contribution(
        &self,
        cfg: &Config,
        two_rows: bool,
        x0_y0_z0_m: Vector3<f64>,
        contribution: &mut SVContribution,
    ) -> Result<VectorContribution, Error> {
        let mut bias_m = 0.0;
        let mut vec = VectorContribution::default();

        let (x0_m, y0_m, z0_m) = (x0_y0_z0_m[0], x0_y0_z0_m[1], x0_y0_z0_m[2]);

        let orbit = self.orbit.ok_or(Error::UnresolvedState)?;
        let pos_vel_m = orbit.to_cartesian_pos_vel() * 1.0E3;

        let (elev_deg, azim_deg) = self.attitude().ok_or(Error::UnresolvedState)?;

        contribution.elevation_deg = elev_deg;
        contribution.azimuth_deg = azim_deg;

        let (sv_x_m, sv_y_m, sv_z_m) = (pos_vel_m[0], pos_vel_m[1], pos_vel_m[2]);

        let mut rho =
            ((sv_x_m - x0_m).powi(2) + (sv_y_m - y0_m).powi(2) + (sv_z_m - z0_m).powi(2)).sqrt();

        rho += self.relativistic_path_range;
        contribution.relativistic_path_range_m = self.relativistic_path_range;

        let (_, range_m) = match cfg.method {
            Method::SPP => {
                let (carrier, pr) = self.best_snr_range_m().ok_or(Error::MissingPseudoRange)?;
                contribution.signal = Signal::Single(carrier);

                (carrier.wavelength(), pr)
            },
            _ => {
                let comb = self
                    .code_if_combination()
                    .ok_or(Error::PseudoRangeCombination)?;

                contribution.signal = Signal::Dual((comb.lhs, comb.rhs));

                let (f1, f2) = (
                    comb.rhs.frequency_hz().powi(2),
                    comb.lhs.frequency_hz().powi(2),
                );

                (SPEED_OF_LIGHT_M_S * (f1 - f2) / f1 / f2, comb.value)
            },
        };

        let cp = match cfg.method {
            Method::PPP => {
                let comb = self
                    .phase_if_combination()
                    .ok_or(Error::PhaseRangeCombination)?;

                Some(comb.value)
            },
            _ => None,
        };

        bias_m -= self.clock_corr.duration.to_seconds() * SPEED_OF_LIGHT_M_S;

        let sys_t = self.system_correction.unwrap_or(Duration::ZERO);

        bias_m += sys_t.to_seconds() * SPEED_OF_LIGHT_M_S;

        if sys_t >= Duration::MIN_POSITIVE {
            debug!(
                "{}({}) - system correction : {}",
                self.epoch, self.sv, sys_t
            );
        }

        bias_m -= self.tgd.to_seconds() * SPEED_OF_LIGHT_M_S;

        if let Some(delay_s) = cfg.externalref_delay_s {
            bias_m -= delay_s * SPEED_OF_LIGHT_M_S;
        }

        for _ in cfg.int_delay.iter() {
            // TODO
        }

        bias_m += self.ionod;

        bias_m += self.tropod;
        contribution.tropo_bias = Some(self.tropod);

        let pr = range_m - rho - bias_m;

        let cp = cp.map(|cp| cp - rho - bias_m);

        if (two_rows || cfg.method == Method::PPP) && cp.is_none() {
            return Err(Error::MissingPhaseRange)?;
        }

        // Measurement weighting: row #1 is code (or phase in single-row PPP),
        // row #2 is phase. Deviations follow the configured elevation model.
        if two_rows {
            vec.row_1 = pr;
            vec.row_2 = cp.unwrap_or_default();
            vec.sigma_1 = cfg.solver.code_sigma_m(Some(elev_deg));
            vec.sigma_2 = cfg.solver.phase_sigma_m(Some(elev_deg));
        } else if cfg.method == Method::PPP {
            vec.row_1 = cp.unwrap_or_default();
            vec.sigma_1 = cfg.solver.phase_sigma_m(Some(elev_deg));
        } else {
            vec.row_1 = pr;
            vec.sigma_1 = cfg.solver.code_sigma_m(Some(elev_deg));
        }

        Ok(vec)
    }

    /// Matrix contribution.
    ///
    /// ## Input
    ///  - i: matrix row
    ///  - cfg: [Config] preset
    ///  - x0_y0_z0: position coordinates as ECEF (m)
    pub(crate) fn ppp_matrix_contribution(
        &self,
        cfg: &Config,
        x0_y0_z0_m: Vector3<f64>,
    ) -> (f64, f64, f64) {
        let (x0_m, y0_m, z0_m) = (x0_y0_z0_m[0], x0_y0_z0_m[1], x0_y0_z0_m[2]);

        let orbit = self.orbit.unwrap_or_else(|| {
            panic!(
                "internal error: {}({}) state not fully resolved!",
                self.epoch, self.sv
            );
        });

        let pos_vel_m = orbit.to_cartesian_pos_vel() * 1.0E3;
        let (sv_x_m, sv_y_m, sv_z_m) = (pos_vel_m[0], pos_vel_m[1], pos_vel_m[2]);

        let mut rho =
            ((x0_m - sv_x_m).powi(2) + (y0_m - sv_y_m).powi(2) + (z0_m - sv_z_m).powi(2)).sqrt();

        if cfg.modeling.relativistic_path_range {
            rho += self.relativistic_path_range;
        }

        let (dx_m, dy_m, dz_m) = (
            (x0_m - sv_x_m) / rho,
            (y0_m - sv_y_m) / rho,
            (z0_m - sv_z_m) / rho,
        );

        (dx_m, dy_m, dz_m)
    }
}

#[cfg(test)]
mod test {
    use crate::{
        prelude::{Config, Epoch, Frame, Method, Orbit},
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
    fn spp_matrix_contribution() {
        let earth_frame = build_earth_frame();
        let t0 = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

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

        let cfg = Config::default().with_navigation_method(Method::SPP);

        let x0_y0_z0_m = Vector3::new(
            ROVER_REFERENCE_COORDS_ECEF_M.0,
            ROVER_REFERENCE_COORDS_ECEF_M.1,
            ROVER_REFERENCE_COORDS_ECEF_M.2,
        );

        let (dx, dy, dz) = rover.ppp_matrix_contribution(&cfg, x0_y0_z0_m);

        let rho = ((ROVER_REFERENCE_COORDS_ECEF_M.0 - e01_position_ecef_m.0).powi(2)
            + (ROVER_REFERENCE_COORDS_ECEF_M.1 - e01_position_ecef_m.1).powi(2)
            + (ROVER_REFERENCE_COORDS_ECEF_M.2 - e01_position_ecef_m.2).powi(2))
        .sqrt();

        let e_i = (
            (ROVER_REFERENCE_COORDS_ECEF_M.0 - e01_position_ecef_m.0) / rho,
            (ROVER_REFERENCE_COORDS_ECEF_M.1 - e01_position_ecef_m.1) / rho,
            (ROVER_REFERENCE_COORDS_ECEF_M.2 - e01_position_ecef_m.2) / rho,
        );

        assert!((dx - e_i.0).abs() < 1E-6, "x error too large");
        assert!((dy - e_i.1).abs() < 1E-6, "y error too large");
        assert!((dz - e_i.2).abs() < 1E-6, "z error too large");
    }

    #[test]
    fn cpp_matrix_contribution() {
        let earth_frame = build_earth_frame();
        let t0 = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

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

        let cfg = Config::default().with_navigation_method(Method::CPP);

        let x0_y0_z0_m = Vector3::new(
            ROVER_REFERENCE_COORDS_ECEF_M.0,
            ROVER_REFERENCE_COORDS_ECEF_M.1,
            ROVER_REFERENCE_COORDS_ECEF_M.2,
        );

        let (dx, dy, dz) = rover.ppp_matrix_contribution(&cfg, x0_y0_z0_m);

        let rho = ((ROVER_REFERENCE_COORDS_ECEF_M.0 - e01_position_ecef_m.0).powi(2)
            + (ROVER_REFERENCE_COORDS_ECEF_M.1 - e01_position_ecef_m.1).powi(2)
            + (ROVER_REFERENCE_COORDS_ECEF_M.2 - e01_position_ecef_m.2).powi(2))
        .sqrt();

        let e_i = (
            (ROVER_REFERENCE_COORDS_ECEF_M.0 - e01_position_ecef_m.0) / rho,
            (ROVER_REFERENCE_COORDS_ECEF_M.1 - e01_position_ecef_m.1) / rho,
            (ROVER_REFERENCE_COORDS_ECEF_M.2 - e01_position_ecef_m.2) / rho,
        );

        assert!((dx - e_i.0).abs() < 1E-6, "x error too large");
        assert!((dy - e_i.1).abs() < 1E-6, "y error too large");
        assert!((dz - e_i.2).abs() < 1E-6, "z error too large");
    }

    #[test]
    fn ppp_matrix_contribution() {
        let earth_frame = build_earth_frame();
        let t0 = Epoch::from_str("2020-06-25T00:00:00 GPST").unwrap();

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

        let cfg = Config::default().with_navigation_method(Method::PPP);

        let x0_y0_z0_m = Vector3::new(
            ROVER_REFERENCE_COORDS_ECEF_M.0,
            ROVER_REFERENCE_COORDS_ECEF_M.1,
            ROVER_REFERENCE_COORDS_ECEF_M.2,
        );

        let (dx, dy, dz) = rover.ppp_matrix_contribution(&cfg, x0_y0_z0_m);

        let rho = ((ROVER_REFERENCE_COORDS_ECEF_M.0 - e01_position_ecef_m.0).powi(2)
            + (ROVER_REFERENCE_COORDS_ECEF_M.1 - e01_position_ecef_m.1).powi(2)
            + (ROVER_REFERENCE_COORDS_ECEF_M.2 - e01_position_ecef_m.2).powi(2))
        .sqrt();

        let e_i = (
            (ROVER_REFERENCE_COORDS_ECEF_M.0 - e01_position_ecef_m.0) / rho,
            (ROVER_REFERENCE_COORDS_ECEF_M.1 - e01_position_ecef_m.1) / rho,
            (ROVER_REFERENCE_COORDS_ECEF_M.2 - e01_position_ecef_m.2) / rho,
        );

        assert!((dx - e_i.0).abs() < 1E-6, "x error too large");
        assert!((dy - e_i.1).abs() < 1E-6, "y error too large");
        assert!((dz - e_i.2).abs() < 1E-6, "z error too large");
    }
}
