use log::error;

use std::collections::HashMap;

use crate::{
    candidate::differences::Differences,
    constants::EARTH_ANGULAR_VEL_RAD,
    prelude::{
        Almanac, Candidate, Config, Duration, EnvironmentalBias, Ephemeris, EphemerisSource, Epoch,
        Frame, Orbit, OrbitSource, Rc, SV, SpacebornBias,
    },
    // smoothing::Smoother,
};

use nalgebra::{Matrix3, Vector3};

pub mod postfit;
pub mod prefit;

pub struct Pool<EPH: EphemerisSource, ORB: OrbitSource, EB: EnvironmentalBias, SB: SpacebornBias> {
    /// [Config]uration preset
    cfg: Config,

    /// Internal [Almanac]
    almanac: Almanac,

    /// ECEF [Frame]
    earth_cef: Frame,

    /// Current [Candidate]s pool
    inner: Vec<Candidate>,

    /// Previous [Candidate]s pool
    past: Vec<Candidate>,

    // /// Measurements [Smoother]
    // smoother: Smoother,
    /// [OrbitSource]
    orb_source: Rc<ORB>,

    /// [EphemerisSource]
    eph_source: Rc<EPH>,

    /// Environmental biases
    env_bias: Rc<EB>,

    /// Spaceborn biases
    space_bias: Rc<SB>,

    /// [Ephemeris] Buffer
    eph_buffer: HashMap<SV, Ephemeris>,

    /// Runtime selected pivot position (ECEF, meters)
    pub pivot_position_ecef_m: Option<(f64, f64, f64)>,

    /// Single [Differences]
    single_differences: Differences,
}

fn orbit_rotation(t: Epoch, dt: Duration, orbit: &Orbit, modeling: bool, frame: Frame) -> Orbit {
    let we = EARTH_ANGULAR_VEL_RAD * dt.to_seconds();
    let (we_sin, we_cos) = we.sin_cos();
    let dcm3 = if modeling {
        Matrix3::new(we_cos, we_sin, 0.0, -we_sin, we_cos, 0.0, 0.0, 0.0, 1.0)
    } else {
        Matrix3::new(1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0)
    };

    let state_m = orbit.to_cartesian_pos_vel() * 1.0E3;
    let position_m = Vector3::new(state_m[0], state_m[1], state_m[2]);
    let position = dcm3 * position_m;

    Orbit::from_position(
        position[0] / 1.0E3,
        position[1] / 1.0E3,
        position[2] / 1.0E3,
        t,
        frame,
    )
}

impl<EPH: EphemerisSource, ORB: OrbitSource, EB: EnvironmentalBias, SB: SpacebornBias>
    Pool<EPH, ORB, EB, SB>
{
    /// Allocate new [Pool]
    pub fn allocate(
        almanac: Almanac,
        cfg: Config,
        earth_cef: Frame,
        eph_source: Rc<EPH>,
        orb_source: Rc<ORB>,
        env_bias: Rc<EB>,
        space_bias: Rc<SB>,
    ) -> Self {
        // let smoother = Smoother::new(cfg.code_smoothing);

        Self {
            cfg,
            earth_cef,
            eph_source,
            orb_source,
            env_bias,
            space_bias,
            // smoother,
            almanac,
            pivot_position_ecef_m: None,
            past: Vec::with_capacity(8),
            inner: Vec::with_capacity(8),
            eph_buffer: HashMap::with_capacity(8),
            single_differences: Default::default(),
        }
    }

    /// Prepare for new epoch
    pub fn new_epoch(&mut self, candidates: &[Candidate]) {
        self.inner = candidates.to_vec();
        self.pivot_position_ecef_m = None;
    }

    /// Returns total number of [Candidate]s
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn candidates(&self) -> &[Candidate] {
        &self.inner
    }

    // pub fn contains(&self, sv: SV) -> bool {
    //     self.inner.iter().filter(|cd| cd.sv == sv).count() > 0
    // }

    pub fn retain<F>(&mut self, f: F)
    where
        F: FnMut(&Candidate) -> bool,
    {
        self.inner.retain(f)
    }

    pub fn retain_mut<F>(&mut self, f: F)
    where
        F: FnMut(&mut Candidate) -> bool,
    {
        self.inner.retain_mut(f)
    }

    /// Determine orbital states
    pub fn orbital_states_fit(&mut self, name: &str) {
        self.inner
            .retain_mut(|cd| match cd.transmission_time(name, &self.cfg) {
                Ok(_) => {
                    // direct state determination
                    let mut determined = false;

                    if let Some(orbit) =
                        &self.orb_source.state_at(cd.tx_epoch, cd.sv, self.earth_cef)
                    {
                        let orbit = orbit_rotation(
                            cd.epoch,
                            cd.signal_time_of_flight(),
                            orbit,
                            self.cfg.modeling.earth_rotation,
                            self.earth_cef,
                        );

                        cd.orbit = Some(orbit);
                        determined = true;
                    }

                    // indirect state determination
                    if !determined {
                        if let Some(eph) = &self.eph_buffer.get(&cd.sv) {
                            if let Some(state) = eph.resolve_state(cd.tx_epoch, self.earth_cef) {
                                let state = orbit_rotation(
                                    cd.epoch,
                                    cd.signal_time_of_flight(),
                                    &state,
                                    self.cfg.modeling.earth_rotation,
                                    self.earth_cef,
                                );

                                cd.orbit = Some(state);
                                determined = true;
                            }
                        }
                    }

                    determined
                },
                Err(e) => {
                    error!("{}({}) {} - tx time error: {}", cd.epoch, cd.sv, name, e);
                    false
                },
            });
    }
}
