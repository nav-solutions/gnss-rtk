use log::{debug, error, info};

use anise::{
    math::Vector3,
    prelude::{Almanac, Frame},
};

use crate::{
    bancroft::Bancroft,
    candidate::Candidate,
    cfg::Config,
    ephemeris::EphemerisSource,
    navigation::{Navigation, apriori::Apriori, solutions::PVTSolution, state::State},
    orbit::OrbitSource,
    pool::Pool,
    prelude::{EnvironmentalBias, Epoch, Error, Rc, SpacebornBias, UserParameters},
    rtk::{NullRTK, RTKBase},
    time::AbsoluteTime,
};

/// [Solver] to resolve [PVTSolution]s.
///
/// ## Generics:
/// - EPH: [EphemerisSource] custom data source. Curreuntly unused: limited to [ORB] only!
/// - ORB: [OrbitSource], custom Orbit data source.
/// - EB: [EnvironmentalBias] implementation using either intenral or custom external model.
/// - SB: [SpacebornBias] external information provider.
/// - TIM: [AbsoluteTime] source for correct absolute time management.
pub struct Solver<
    EPH: EphemerisSource,
    ORB: OrbitSource,
    EB: EnvironmentalBias,
    SB: SpacebornBias,
    TIM: AbsoluteTime,
> {
    /// Solver [Config]uration
    pub cfg: Config,

    /// [Frame]
    earth_cef: Frame,

    /// Rover pool
    rover_pool: Pool<EPH, ORB, EB, SB>,

    /// Base pool
    base_pool: Pool<EPH, ORB, EB, SB>,

    /// PPP prefit

    /// [Navigation] solver
    navigation: Navigation,

    /// [AbsoluteTime] implementation
    absolute_time: TIM,

    /// Possible initial position
    initial_ecef_m: Option<Vector3>,
}

impl<
    EPH: EphemerisSource,
    ORB: OrbitSource,
    EB: EnvironmentalBias,
    SB: SpacebornBias,
    TIM: AbsoluteTime,
> Solver<EPH, ORB, EB, SB, TIM>
{
    /// Creates a new [Solver] for either direct or differential navigation,
    /// with possible apriori knowledge.
    ///
    /// ## Input
    /// - almanac: provided valid [Almanac]
    /// - earth_cef: [Frame] that must be an ECEF
    /// - cfg: solver [Config]uration
    /// - eph_source: custom [EphemerisSource], unavailable right now,
    /// tie to null.
    /// - orbit_source: custom [OrbitSource] implementation,
    /// wrapped in a Rc<RefCell<>> which allows the solver
    /// and the orbital provider to live in the same thread.
    /// - absolute_time: external [AbsoluteTime] implementation.
    /// - bias: [Bias] model implementation
    /// - state_ecef_m: provide initial state as ECEF 3D coordinates,
    /// otherwise we will have to figure them.
    pub fn new(
        almanac: Almanac,
        earth_cef: Frame,
        cfg: Config,
        eph_source: Rc<EPH>,
        orbit_source: Rc<ORB>,
        spaceborn_biases: Rc<SB>,
        environmental_biases: Rc<EB>,
        absolute_time: TIM,
        state_ecef_m: Option<(f64, f64, f64)>,
    ) -> Self {
        let initial_ecef_m = match state_ecef_m {
            Some((x0, y0, z0)) => Some(Vector3::new(x0, y0, z0)),
            _ => None,
        };

        let navigation = Navigation::new(&cfg, earth_cef);

        let rover_pool = Pool::allocate(
            almanac.clone(),
            cfg.clone(),
            earth_cef,
            eph_source.clone(),
            orbit_source.clone(),
            environmental_biases.clone(),
            spaceborn_biases.clone(),
        );

        let base_pool = Pool::allocate(
            almanac,
            cfg.clone(),
            earth_cef,
            eph_source,
            orbit_source,
            environmental_biases,
            spaceborn_biases,
        );

        Self {
            earth_cef,
            navigation,
            rover_pool,
            base_pool,
            absolute_time,
            initial_ecef_m,
            cfg: cfg.clone(),
        }
    }

    /// Creates a new [Solver] to operate without a priori knowledge
    /// and perform a complete survey from scratch. Refer to [Self::new]
    /// for more information.
    pub fn new_survey(
        almanac: Almanac,
        earth_cef: Frame,
        cfg: Config,
        eph_source: Rc<EPH>,
        orbit_source: Rc<ORB>,
        spaceborn_biases: Rc<SB>,
        environmental_biases: Rc<EB>,
        absolute_time: TIM,
    ) -> Self {
        Self::new(
            almanac,
            earth_cef,
            cfg,
            eph_source,
            orbit_source,
            spaceborn_biases,
            environmental_biases,
            absolute_time,
            None,
        )
    }

    /// [PVTSolution] solving attempt using PPP technique (no reference).
    /// Unlike RTK resolution attempt, PPP solving will resolve both
    /// the position and the local clock state, but it is
    /// a less accurate and much more complex process.
    ///
    /// ## Input
    /// - epoch: [Epoch] of measurement
    /// - params: [UserParameters]
    /// - candidates: proposed [Candidate]s (= measurements)
    /// - rtk_base: possible [RTKBase] we will connect to
    ///
    /// ## Output
    /// - success: [PVTSolution]
    /// - failure: [Error]
    pub fn ppp(
        &mut self,
        epoch: Epoch,
        params: UserParameters,
        candidates: &[Candidate],
    ) -> Result<PVTSolution, Error> {
        let null_rtk = NullRTK {};
        let solution = self.solve(epoch, params, candidates, &null_rtk, false)?;
        Ok(solution)
    }

    /// [PVTSolution] solving attempt using RTK technique and a single remote
    /// reference site that implements [RTKBase].
    /// You may catch RTK solving issues (for example, network loss)
    /// and retry using [Self::ppp], because this object is compatible
    /// with both. For this function to converge, the remote
    /// site and rover proposals must agree in their content (enough matches)
    /// and signals must fulfill the navigation [Method] being used.
    ///
    /// NB: unlike PPP solving, RTK solving will only resolve the spatial
    /// state, the clock state can only remain undetermined. If you are
    /// interested by both, you will either have to restrict to [Self::ppp_run],
    /// or do a PPP (time only) run after the RTK run.
    ///
    /// ## Input
    /// - epoch: [Epoch] of measurement from the rover reported by the rover.
    /// The remote site must match the sampling instant closely, reducing
    /// the time-difference (RTK aging).
    /// - profile: rover [UserProfile]
    /// - candidates: rover measurements, wrapped as [Candidate]s.
    /// - rtk_base: remote reference site that implements [RTKBase]reference.
    ///
    /// ## Output
    /// - [PVTSolution].
    pub fn rtk<RTK: RTKBase>(
        &mut self,
        epoch: Epoch,
        params: UserParameters,
        candidates: &[Candidate],
        rtk_base: &RTK,
    ) -> Result<PVTSolution, Error> {
        self.solve(epoch, params, candidates, rtk_base, true)
    }

    /// [PVTSolution] solving attempt.
    fn solve<RTK: RTKBase>(
        &mut self,
        epoch: Epoch,
        params: UserParameters,
        pool: &[Candidate],
        rtk_base: &RTK,
        uses_rtk: bool,
    ) -> Result<PVTSolution, Error> {
        let rtk_base_name = rtk_base.name();
        let min_required = self.min_sv_required(uses_rtk);

        if pool.len() < min_required {
            // no need to proceed further
            return Err(Error::NotEnoughCandidates);
        }

        self.rover_pool.new_epoch(pool);

        if uses_rtk {
            let observations = rtk_base.observe(epoch);
            info!("{epoch} - using remote {rtk_base_name} reference");
            self.base_pool.new_epoch(&observations);
        }

        self.rover_pool.pre_fit("rover", &self.absolute_time);

        if uses_rtk {
            self.base_pool.pre_fit(&rtk_base_name, &self.absolute_time);
        }

        if self.rover_pool.len() < min_required {
            return Err(Error::NotEnoughPreFitCandidates);
        }

        let epoch = if epoch.time_scale == self.cfg.timescale {
            epoch
        } else {
            let corrected = self
                .absolute_time
                .epoch_correction(epoch, self.cfg.timescale);

            debug!(
                "{} - |{}-{}| corrected sampling: {}",
                epoch, epoch.time_scale, self.cfg.timescale, corrected
            );

            corrected
        };

        self.rover_pool.orbital_states_fit("rover");

        if uses_rtk {
            self.base_pool.orbital_states_fit(&rtk_base_name);
        }

        // current state
        let state = if self.navigation.is_initialized() {
            self.navigation.state.with_epoch(epoch)
        } else {
            match self.initial_ecef_m {
                Some(x0_y0_z0_m) => {
                    let apriori = Apriori::from_ecef_m(x0_y0_z0_m, epoch, self.earth_cef);

                    let state = State::from_apriori(&apriori).unwrap_or_else(|e| {
                        panic!("Solver initial preset is physically incorrect: {e}");
                    });

                    debug!("{epoch} - initial state: {state}");
                    state
                },
                None => {
                    let solver = Bancroft::new(self.rover_pool.candidates())?;
                    let solution = solver.resolve()?;
                    let x0_y0_z0_m = Vector3::new(solution[0], solution[1], solution[2]);

                    let apriori = Apriori::from_ecef_m(x0_y0_z0_m, epoch, self.earth_cef);

                    let state = State::from_apriori(&apriori).unwrap_or_else(|e| {
                        panic!("Solver failed to initialize itself. Physical error: {e}");
                    });

                    debug!("{epoch} - initial state: {state}");
                    state
                },
            }
        };

        // rover post-fit
        self.rover_pool.post_fit("rover", &state).map_err(|e| {
            error!("{epoch} rover postfit error {e}");
            Error::PostfitPrenav
        })?;

        let double_differences = if uses_rtk {
            // base post-fit
            self.base_pool
                .post_fit(&rtk_base_name, &state)
                .map_err(|e| {
                    error!("{epoch} {rtk_base_name} postfit error: {e}");
                    Error::PostfitPrenav
                })?;

            // special RTK post-fit
            match self.rover_pool.rtk_post_fit(&mut self.base_pool) {
                Ok(double_differences) => Some(double_differences),
                Err(e) => {
                    error!("{epoch} - rtk post-fit error: {e}");
                    return Err(e);
                },
            }
        } else {
            None
        };

        let pool_size = self.rover_pool.len();

        if pool_size < min_required {
            return Err(Error::NotEnoughPostFitCandidates);
        }

        // Solving attempt
        match self.navigation.solve(
            epoch,
            params,
            &state,
            self.rover_pool.candidates(),
            pool_size,
            uses_rtk,
            rtk_base,
            &self.rover_pool.pivot_position_ecef_m,
            &double_differences,
        ) {
            Ok(_) => {
                info!("{epoch} - sucess");
            },
            Err(e) => {
                error!("{epoch} - iteration failure: {e}");
                return Err(e);
            },
        }

        // Form P.V.T solution
        let solution = PVTSolution::new(
            epoch,
            uses_rtk,
            &self.navigation.state,
            &self.navigation.dop,
            &self.navigation.sv,
        );

        // Special "open loop" option
        if self.cfg.solver.open_loop {
            self.navigation.state = state;
        }

        Ok(solution)
    }

    /// Reset this [Solver].
    pub fn reset(&mut self) {
        self.navigation.reset();
    }

    /// Returns minimal requirement for current preset
    fn min_sv_required(&self, uses_rtk: bool) -> usize {
        let mut min_sv = 4;

        if self.navigation.is_initialized() && self.cfg.fixed_altitude.is_some() {
            min_sv -= 1;
        }

        if uses_rtk {
            min_sv += 1;
        }

        min_sv
    }
}
