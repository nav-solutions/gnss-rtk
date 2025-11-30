use nalgebra::{DVector, DimName, U3, U4};

use anise::{
    astro::PhysicsResult,
    math::{Vector3, Vector6},
    prelude::{Epoch, Frame},
};

use crate::{
    constants::SPEED_OF_LIGHT_M_S,
    navigation::{Apriori, Navigation},
    prelude::Orbit,
};

#[derive(Debug, Clone, PartialEq)]
pub struct State {
    /// [Epoch] of resolution
    pub epoch: Epoch,

    /// [DVector]
    pub x: DVector<f64>,

    /// x_amb [DVector]
    pub x_amb: DVector<f64>,

    /// Clock drift (s.s⁻¹)
    pub clock_drift_s_s: f64,

    /// Geodeticy position (ddeg, ddeg, km above mean sea level)
    pub lat_long_alt_deg_deg_km: (f64, f64, f64),
}

impl Default for State {
    /// Builds a default Null [State].
    fn default() -> Self {
        Self {
            epoch: Default::default(),
            x_amb: Default::default(),
            x: DVector::<f64>::zeros(U4::DIM),
            clock_drift_s_s: Default::default(),
            lat_long_alt_deg_deg_km: Default::default(),
        }
    }
}

impl std::fmt::Display for State {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let position_ecef_m = self.to_position_ecef_m();
        let (offset, drift) = self.clock_profile_s();

        writeln!(
            f,
            "{position_ecef_m} dt={offset:.11E}s drift={drift:.11E}s/s"
        )?;

        if self.x_amb.nrows() > 0 {
            writeln!(f, "ambiguities={}", self.x_amb)?;
        }

        Ok(())
    }
}

impl State {
    /// Create new [State] initialized from [Apriori]
    pub fn from_apriori(apriori: &Apriori) -> PhysicsResult<Self> {
        let orbit = apriori.to_orbit();
        Self::from_orbit(&orbit)
    }

    /// Create a new [State] with updated [Epoch]
    pub fn with_epoch(&self, epoch: Epoch) -> Self {
        let mut s = self.clone();
        s.epoch = epoch;
        s
    }

    pub fn resize_mut(&mut self, ndf: usize) {
        self.x.resize_vertically_mut(ndf, 0.0);
    }

    pub fn resize_ambiguities_mut(&mut self, ndf: usize) {
        self.x_amb.resize_vertically_mut(ndf, 0.0);
    }

    /// Create new [State] from [Orbit]al solution.
    pub fn from_orbit(orbit: &Orbit) -> PhysicsResult<Self> {
        let pos_vel_m = orbit.to_cartesian_pos_vel() * 1.0E3;
        let latlongalt = orbit.latlongalt()?;

        let mut x = DVector::<f64>::zeros(U4::DIM);

        for i in 0..U3::DIM {
            x[i] = pos_vel_m[i];
        }

        Ok(Self {
            x,
            epoch: orbit.epoch,
            clock_drift_s_s: 0.0_f64,
            x_amb: Default::default(),
            lat_long_alt_deg_deg_km: latlongalt,
        })
    }

    /// Returns position in ECEF meters as [Vector3]
    pub fn to_position_ecef_m(&self) -> Vector3 {
        Vector3::new(self.x[0], self.x[1], self.x[2])
    }

    /// Returns position and velocity in ECEF meters as [Vector6]
    pub fn to_position_velocity_ecef_m(&self) -> Vector6 {
        Vector6::new(self.x[0], self.x[1], self.x[2], 0.0, 0.0, 0.0)
    }

    /// Returns estimated clock (offset, drift) in seconds and s.s⁻¹.
    pub fn clock_profile_s(&self) -> (f64, f64) {
        if self.x.nrows() > Navigation::clock_index() {
            (self.x[Navigation::clock_index()], self.clock_drift_s_s)
        } else {
            (0.0_f64, 0.0_f64)
        }
    }

    /// Converts [State] to [Orbit]
    pub fn to_orbit(&self, frame: Frame) -> Orbit {
        let pos_vel_km_s = self.to_position_velocity_ecef_m() / 1.0E3;
        Orbit::from_cartesian_pos_vel(pos_vel_km_s, self.epoch, frame)
    }

    /// Spatial [State] correction with mutable access.
    pub fn spatial_correction_mut(
        &mut self,
        frame: Frame,
        dx: (f64, f64, f64),
    ) -> PhysicsResult<()> {
        self.x[0] += dx.0;
        self.x[1] += dx.1;
        self.x[2] += dx.2;

        // update attitude
        let new_orbit = self.to_orbit(frame);
        self.lat_long_alt_deg_deg_km = new_orbit.latlongalt()?;

        Ok(())
    }

    /// Temporal [State] correction, with mutable access.
    pub fn temporal_correction_mut(&mut self, dt: f64) {
        // let dt = (pending_t - self.epoch).to_seconds();
        // if dt > 0.0 && ndf > clock_index {
        //     self.clock_drift_s_s =
        //         ((dx[clock_index] / SPEED_OF_LIGHT_M_S) - self.x[clock_index]) / dt;
        // }

        self.x[Navigation::clock_index()] = dt / SPEED_OF_LIGHT_M_S;
    }

    /// Temporal update
    pub fn postfit_update_mut(&mut self, frame: Frame, dx: &DVector<f64>) -> PhysicsResult<()> {
        // TODO: velocity
        for i in 0..U3::DIM {
            self.x[i] = dx[i];
        }

        let new_orbit = self.to_orbit(frame);
        self.lat_long_alt_deg_deg_km = new_orbit.latlongalt()?;

        Ok(())
    }

    // /// Compute residual between this [State] and other [State]
    // pub(crate) fn residual(&self, rhs: &State<D>) -> OVector<f64, D> {
    //     self.x - rhs.x
    // }
}

#[cfg(test)]
mod test {

    use crate::{
        navigation::{Apriori, State},
        prelude::{Duration, Frame, Orbit, SPEED_OF_LIGHT_M_S},
        tests::ROVER_REFERENCE_COORDS_ECEF_M,
    };

    use anise::math::Vector3;

    // use nalgebra::{DimName, U4};

    use rstest::*;

    #[fixture]
    fn build_earth_frame() -> Frame {
        use crate::tests::earth_frame;
        earth_frame()
    }

    #[fixture]
    fn build_reference_apriori() -> Apriori {
        use crate::tests::rover_reference_apriori_at_ref_epoch;
        rover_reference_apriori_at_ref_epoch()
    }

    #[fixture]
    fn build_reference_orbit() -> Orbit {
        use crate::tests::{earth_frame, rover_reference_orbit_at_ref_epoch};
        let earth_frame = earth_frame();
        rover_reference_orbit_at_ref_epoch(earth_frame)
    }

    #[test]
    fn state_u4() {
        let earth_frame = build_earth_frame();
        let apriori = build_reference_apriori();
        let reference_orbit = build_reference_orbit();

        let initial_state = State::from_apriori(&apriori).unwrap_or_else(|e| {
            panic!("Failed to build initial state from reference apriori: {e}")
        });

        // verify that orbital construction works correctly
        let from_orbit = State::from_orbit(&reference_orbit)
            .unwrap_or_else(|e| panic!("Failed to build initial state from reference orbit: {e}"));

        assert_eq!(initial_state, from_orbit);

        // verify initial state
        let position_ecef_m = initial_state.to_position_ecef_m();

        assert_eq!(
            position_ecef_m,
            Vector3::new(
                ROVER_REFERENCE_COORDS_ECEF_M.0,
                ROVER_REFERENCE_COORDS_ECEF_M.1,
                ROVER_REFERENCE_COORDS_ECEF_M.2
            )
        );

        let position_velocity_m_s = initial_state.to_position_velocity_ecef_m();

        assert_eq!(
            (
                position_velocity_m_s[0],
                position_velocity_m_s[1],
                position_velocity_m_s[2]
            ),
            (
                ROVER_REFERENCE_COORDS_ECEF_M.0,
                ROVER_REFERENCE_COORDS_ECEF_M.1,
                ROVER_REFERENCE_COORDS_ECEF_M.2
            ),
            "initial position error!"
        );

        assert_eq!(
            (
                position_velocity_m_s[3],
                position_velocity_m_s[4],
                position_velocity_m_s[5],
            ),
            (0.0, 0.0, 0.0),
            "dynamics not modeled in yet!",
        );

        assert_eq!(
            initial_state.clock_profile_s(),
            (0.0, 0.0),
            "clock state should not be initialized!"
        );

        // verify +null is null
        let mut state = initial_state.clone();
        let null_dx = (0.0, 0.0, 0.0);

        state
            .spatial_correction_mut(earth_frame, null_dx)
            .unwrap_or_else(|e| panic!("Failed to apply null state correction! {e}"));

        assert_eq!(
            initial_state, state,
            "null correction should not modify internal state!"
        );

        // Test correction (2)
        let mut state = initial_state.clone();

        let (dx, dy, dz) = (1.0, 2.0, 3.0);

        state
            .spatial_correction_mut(earth_frame, (dx, dy, dz))
            .unwrap_or_else(|e| panic!("Failed to apply state correction! {e}"));

        assert_eq!(
            state.epoch, initial_state.epoch,
            "Epoch should have been preserved"
        );

        state.temporal_correction_mut(4.0);

        let position_ecef_m = state.to_position_ecef_m();

        assert_eq!(
            position_ecef_m,
            Vector3::new(
                ROVER_REFERENCE_COORDS_ECEF_M.0 + 1.0,
                ROVER_REFERENCE_COORDS_ECEF_M.1 + 2.0,
                ROVER_REFERENCE_COORDS_ECEF_M.2 + 3.0,
            ),
            "invalid spatial state",
        );

        let (clock_offset_s, clock_drift_s) = state.clock_profile_s();

        assert_eq!(
            clock_offset_s,
            4.0 / SPEED_OF_LIGHT_M_S,
            "invalid temporal state"
        );
        assert_eq!(clock_drift_s, 0.0, "clock drift should have been preserved");

        // Test correction (2)
        let mut state = initial_state.clone();
        let new_t = state.epoch + Duration::from_seconds(30.0);

        let (dx, dy, dz) = (1.0, 2.0, 3.0);

        state
            .spatial_correction_mut(earth_frame, (dx, dy, dz))
            .unwrap_or_else(|e| panic!("Failed to apply null state correction! {e}"));

        state.temporal_correction_mut(4.0);

        // assert!(
        //     state.epoch > initial_state.epoch,
        //     "Epoch should have been updated!"
        // );

        // assert_eq!(
        //     (state.epoch - initial_state.epoch),
        //     Duration::from_seconds(30.0),
        //     "Invalid sampling update"
        // );

        let position_ecef_m = state.to_position_ecef_m();

        assert_eq!(
            position_ecef_m,
            Vector3::new(
                ROVER_REFERENCE_COORDS_ECEF_M.0 + 1.0,
                ROVER_REFERENCE_COORDS_ECEF_M.1 + 2.0,
                ROVER_REFERENCE_COORDS_ECEF_M.2 + 3.0,
            ),
            "invalid spatial state correction",
        );

        let (clock_offset_s, clock_drift_s) = state.clock_profile_s();

        assert_eq!(clock_offset_s, 4.0 / SPEED_OF_LIGHT_M_S);

        // assert_eq!(
        //     clock_drift_s,
        //     4.0 / SPEED_OF_LIGHT_M_S / Duration::from_seconds(30.0).to_seconds(),
        //     "invalid clock drift update!"
        // );
    }
}
