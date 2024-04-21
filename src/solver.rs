//! PVT solver
use std::collections::HashMap;

use gnss::prelude::SV;
use hifitime::Unit;
use log::{debug, error, warn};
use map_3d::deg2rad;
use nalgebra::{Matrix3, Vector3};
use thiserror::Error;

use nyx::{
    cosmic::{
        eclipse::{eclipse_state, EclipseState},
        Orbit, SPEED_OF_LIGHT,
    },
    md::prelude::{Arc, Cosm, Frame},
};

use crate::{
    apriori::AprioriPosition,
    bias::{IonosphereBias, TroposphereBias},
    candidate::Candidate,
    cfg::{Config, Method},
    navigation::{
        solutions::validator::Validator as SolutionValidator, Input as NavigationInput, Navigation,
        PVTSolution, PVTSolutionType,
    },
    prelude::{Duration, Epoch},
};

#[derive(Debug, Clone, Error)]
pub enum Error {
    #[error("need more candidates to resolve a {0} a solution")]
    NotEnoughInputCandidates(PVTSolutionType),
    #[error("not enough candidates fit criteria")]
    NotEnoughFittingCandidates,
    #[error("failed to form matrix")]
    MatrixError,
    #[error("failed to invert matrix")]
    MatrixInversionError,
    #[error("reolved NaN: invalid input matrix")]
    TimeIsNan,
    #[error("undefined apriori position")]
    UndefinedAprioriPosition,
    #[error("missing pseudo range observation")]
    MissingPseudoRange,
    #[error("cannot form signal combination: missing dual freq signals")]
    PseudoRangeCombination,
    #[error("at least one pseudo range observation is mandatory")]
    NeedsAtLeastOnePseudoRange,
    #[error("failed to model or measure ionospheric delay")]
    MissingIonosphericDelayValue,
    #[error("unresolved state: interpolation should have passed")]
    UnresolvedState,
    #[error("unable to form signal combination")]
    SignalRecombination,
    #[error("physical non sense: rx prior tx")]
    PhysicalNonSenseRxPriorTx,
    #[error("physical non sense: t_rx is too late")]
    PhysicalNonSenseRxTooLate,
    #[error("invalidated solution")]
    InvalidatedSolution,
    // Kalman filter bad op: should never happen
    #[error("uninitialized kalman filter!")]
    UninitializedKalmanFilter,
}

/// Interpolation result (state vector) that needs to be
/// resolved for every single candidate.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub struct InterpolationResult {
    /// Elevation compared to reference position and horizon in [°]
    pub elevation: f64,
    /// Azimuth compared to reference position and magnetic North in [°]
    pub azimuth: f64,
    /// APC Position vector in [m] ECEF
    pub position: Vector3<f64>,
    // Velocity vector in [m/s] ECEF that we calculated ourselves
    velocity: Option<Vector3<f64>>,
}

impl InterpolationResult {
    /// Builds InterpolationResults from an Antenna Phase Center (APC) position
    /// as ECEF [m] coordinates
    pub fn from_apc_position(position: (f64, f64, f64)) -> Self {
        Self {
            velocity: None,
            azimuth: 0.0_f64,
            elevation: 0.0_f64,
            position: Vector3::<f64>::new(position.0, position.1, position.2),
        }
    }
    /// Builds Self with given SV (elevation, azimuth) attitude
    pub fn with_elevation_azimuth(&self, elev_azim: (f64, f64)) -> Self {
        Self {
            azimuth: elev_azim.1,
            elevation: elev_azim.0,
            position: self.position,
            velocity: self.velocity,
        }
    }
    pub(crate) fn velocity(&self) -> Option<Vector3<f64>> {
        self.velocity
    }
    pub(crate) fn orbit(&self, dt: Epoch, frame: Frame) -> Orbit {
        let p = self.position;
        let v = self.velocity().unwrap_or_default();
        Orbit::cartesian(
            p[0] / 1000.0,
            p[1] / 1000.0,
            p[2] / 1000.0,
            v[0] / 1000.0,
            v[1] / 1000.0,
            v[2] / 1000.0,
            dt,
            frame,
        )
    }
}

/// PVT Solver.
/// I: Interpolated SV APC coordinates interface.
/// You are required to provide APC coordinates at requested ("t", "sv"),
/// expressed in meters [ECEF], for this to proceed.
#[derive(Debug, Clone)]
pub struct Solver<I>
where
    I: Fn(Epoch, SV, usize) -> Option<InterpolationResult>,
{
    /// Solver parametrization
    pub cfg: Config,
    /// apriori position
    pub apriori: AprioriPosition,
    /// Interpolated SV state.
    pub interpolator: I,
    /* Cosmic model */
    cosmic: Arc<Cosm>,
    /*
     * Solid / Earth body frame.
     */
    earth_frame: Frame,
    /*
     * Sun / Star body frame
     */
    sun_frame: Frame,
    // navigator
    nav: Navigation,
    /* prev. solution for internal logic */
    prev_pvt: Option<(Epoch, PVTSolution)>,
    /* prev. state vector for internal velocity determination */
    prev_sv_state: HashMap<SV, (Epoch, Vector3<f64>)>,
}

impl<I: std::ops::Fn(Epoch, SV, usize) -> Option<InterpolationResult>> Solver<I> {
    pub fn new(cfg: &Config, apriori: AprioriPosition, interpolator: I) -> Result<Self, Error> {
        let cosmic = Cosm::de438();
        let sun_frame = cosmic.frame("Sun J2000");
        let earth_frame = cosmic.frame("EME2000");
        /*
         * print some infos on latched config
         */
        if cfg.method == Method::SPP && cfg.min_sv_sunlight_rate.is_some() {
            warn!("Eclipse filter is not meaningful in SPP mode");
        }
        if cfg.modeling.relativistic_path_range {
            warn!("Relativistic path range cannot be modeled at the moment");
        }
        Ok(Self {
            cosmic,
            sun_frame,
            earth_frame,
            apriori,
            interpolator,
            cfg: cfg.clone(),
            prev_sv_state: HashMap::new(),
            nav: Navigation::new(cfg.solver.filter),
            prev_pvt: Option::<(Epoch, PVTSolution)>::None,
        })
    }
    /// Try to resolve a PVTSolution at desired "t".
    /// "t": sampling instant.
    /// "solution": desired PVTSolutionType.
    /// "pool": List of candidates.
    /// iono_bias: possible IonosphereBias if you can provide such info.
    /// tropo_bias: possible TroposphereBias if you can provide such info.
    pub fn resolve(
        &mut self,
        t: Epoch,
        pvt_type: PVTSolutionType,
        pool: Vec<Candidate>,
        iono_bias: &IonosphereBias,
        tropo_bias: &TroposphereBias,
    ) -> Result<(Epoch, PVTSolution), Error> {
        let min_required = Self::min_required(pvt_type, &self.cfg);

        if pool.len() < min_required {
            return Err(Error::NotEnoughInputCandidates(pvt_type));
        }

        let (x0, y0, z0) = (
            self.apriori.ecef.x,
            self.apriori.ecef.y,
            self.apriori.ecef.z,
        );

        let (lat_ddeg, lon_ddeg, altitude_above_sea_m) = (
            self.apriori.geodetic.x,
            self.apriori.geodetic.y,
            self.apriori.geodetic.z,
        );

        let method = self.cfg.method;
        let modeling = self.cfg.modeling;
        let solver_opts = &self.cfg.solver;
        let _filter = solver_opts.filter;
        let interp_order = self.cfg.interp_order;

        /* apply signal quality and condition filters */
        let pool: Vec<Candidate> = pool
            .into_iter()
            .filter_map(|cd| match method {
                Method::SPP => {
                    let pr = cd.prefered_pseudorange()?;
                    if let Some(min_snr) = self.cfg.min_snr {
                        let snr = pr.snr?;
                        if snr < min_snr {
                            None
                        } else {
                            Some(cd)
                        }
                    } else {
                        Some(cd)
                    }
                },
                Method::CodePPP => {
                    if cd.code_ppp_compatible() {
                        Some(cd)
                    } else {
                        debug!("{:?} ({}) missing secondary frequency", cd.t, cd.sv);
                        None
                    }
                },
            })
            .collect();

        /* interpolate positions */
        let mut pool: Vec<Candidate> = pool
            .iter()
            .filter_map(|cd| match cd.transmission_time(&self.cfg) {
                Ok((t_tx, dt_tx)) => {
                    debug!("{:?} ({}) : signal propagation {}", cd.t, cd.sv, dt_tx);
                    let interpolated = (self.interpolator)(t_tx, cd.sv, interp_order)?;

                    let min_elev = self.cfg.min_sv_elev.unwrap_or(0.0_f64);
                    let min_azim = self.cfg.min_sv_azim.unwrap_or(0.0_f64);
                    let max_azim = self.cfg.max_sv_azim.unwrap_or(360.0_f64);

                    if interpolated.elevation < min_elev {
                        debug!(
                            "{:?} ({}) - {:?} rejected : below elevation mask",
                            cd.t, cd.sv, interpolated
                        );
                        None
                    } else if interpolated.azimuth < min_azim {
                        debug!(
                            "{:?} ({}) - {:?} rejected : below azimuth mask",
                            cd.t, cd.sv, interpolated
                        );
                        None
                    } else if interpolated.azimuth > max_azim {
                        debug!(
                            "{:?} ({}) - {:?} rejected : above azimuth mask",
                            cd.t, cd.sv, interpolated
                        );
                        None
                    } else {
                        let mut cd = cd.clone();
                        let interpolated =
                            Self::rotate_position(modeling.earth_rotation, interpolated, dt_tx);
                        let interpolated = self.velocities(t_tx, cd.sv, interpolated);
                        cd.t_tx = t_tx;
                        debug!("{:?} ({}) : {:?}", cd.t, cd.sv, interpolated);
                        cd.state = Some(interpolated);
                        Some(cd)
                    }
                },
                Err(e) => {
                    error!("{} - transmision time error: {:?}", cd.sv, e);
                    None
                },
            })
            .collect();

        /*
         * Update internal state
         */
        for cd in pool.iter_mut() {
            if modeling.relativistic_clock_bias {
                /*
                 * following calculations need inst. velocity
                 */
                let state = cd.state.unwrap();
                if state.velocity.is_some() {
                    const EARTH_SEMI_MAJOR_AXIS_WGS84: f64 = 6378137.0_f64;
                    const EARTH_GRAVITATIONAL_CONST: f64 = 3986004.418 * 10.0E8;
                    let orbit = state.orbit(cd.t_tx, self.earth_frame);
                    let ea_rad = deg2rad(orbit.ea_deg());
                    let gm = (EARTH_SEMI_MAJOR_AXIS_WGS84 * EARTH_GRAVITATIONAL_CONST).sqrt();
                    let bias = -2.0_f64 * orbit.ecc() * ea_rad.sin() * gm
                        / SPEED_OF_LIGHT
                        / SPEED_OF_LIGHT
                        * Unit::Second;
                    debug!("{:?} ({}) : relativistic clock bias: {}", cd.t, cd.sv, bias);
                    cd.clock_corr += bias;
                }
            }

            self.prev_sv_state
                .insert(cd.sv, (cd.t_tx, cd.state.unwrap().position));
        }

        /* apply eclipse filter (if need be) */
        if let Some(min_rate) = self.cfg.min_sv_sunlight_rate {
            let mut nb_removed: usize = 0;
            for idx in 0..pool.len() {
                let state = pool[idx - nb_removed].state.unwrap(); // infaillible
                let orbit = state.orbit(pool[idx - nb_removed].t, self.earth_frame);
                let state = eclipse_state(&orbit, self.sun_frame, self.earth_frame, &self.cosmic);
                let eclipsed = match state {
                    EclipseState::Umbra => true,
                    EclipseState::Visibilis => false,
                    EclipseState::Penumbra(r) => r < min_rate,
                };
                if eclipsed {
                    debug!(
                        "{:?} ({}): dropped - eclipsed by Earth",
                        pool[idx - nb_removed].t,
                        pool[idx - nb_removed].sv
                    );
                    let _ = pool.swap_remove(idx - nb_removed);
                    nb_removed += 1;
                }
            }
        }

        // adapt to navigation
        pool.sort_by(|cd_a, cd_b| {
            let state_a = cd_a.state.unwrap();
            let state_b = cd_b.state.unwrap();
            state_b.elevation.partial_cmp(&state_a.elevation).unwrap()
        });

        match pvt_type {
            PVTSolutionType::TimeOnly => {
                let mut index = 0;
                pool.retain(|_| {
                    index += 1;
                    index == 1
                });
            },
            _ => {
                let mut index = 0;
                pool.retain(|_| {
                    index += 1;
                    index < min_required + 1
                });
            },
        }

        if pool.len() != min_required {
            return Err(Error::NotEnoughFittingCandidates);
        }

        let input = match NavigationInput::new(
            (x0, y0, z0),
            (lat_ddeg, lon_ddeg, altitude_above_sea_m),
            &self.cfg,
            &pool,
            iono_bias,
            tropo_bias,
        ) {
            Ok(input) => input,
            Err(e) => {
                error!("Failed to form navigation matrix: {}", e);
                return Err(Error::MatrixError);
            },
        };

        let output = match self.nav.resolve(&input) {
            Ok(output) => output,
            Err(e) => {
                error!("Failed to resolve: {}", e);
                return Err(Error::MatrixError);
            },
        };

        let validator = SolutionValidator::new(&self.apriori.ecef, &pool, &input, &output);

        match validator.validate(solver_opts) {
            Ok(_) => self.nav.validate(),
            Err(e) => {
                error!("solution invalidated - {}", e);
                return Err(Error::InvalidatedSolution);
            },
        };

        let x = output.state.estimate();

        let mut solution = PVTSolution {
            q: output.q,
            gdop: output.gdop,
            tdop: output.tdop,
            pdop: output.pdop,
            vel: Vector3::<f64>::default(),
            pos: Vector3::new(x[0], x[1], x[2]),
            sv: input.sv.clone(),
            dt: x[3] / SPEED_OF_LIGHT,
        };

        let _to_discard = false;

        if let Some((prev_t, prev_pvt)) = &self.prev_pvt {
            solution.vel = (solution.pos - prev_pvt.pos) / (t - *prev_t).to_seconds();
        }

        if self.prev_pvt.is_none() {
            // always discard 1st solution
            self.prev_pvt = Some((t, solution.clone()));
            return Err(Error::InvalidatedSolution);
        }

        self.prev_pvt = Some((t, solution.clone()));

        Self::rework_solution(&mut solution, pvt_type, &self.apriori, &self.cfg);
        Ok((t, solution))
    }
    /*
     * Returns nb of vehicles we need to gather
     */
    fn min_required(solution: PVTSolutionType, cfg: &Config) -> usize {
        match solution {
            PVTSolutionType::TimeOnly => 1,
            _ => {
                let mut n = 4;
                if cfg.fixed_altitude.is_some() {
                    n -= 1;
                }
                n
            },
        }
    }
    /*
     * Apply appropriate adjustments
     */
    fn rotate_position(
        rotate: bool,
        interpolated: InterpolationResult,
        dt_tx: Duration,
    ) -> InterpolationResult {
        let mut reworked = interpolated;
        let rot = if rotate {
            const EARTH_OMEGA_E_WGS84: f64 = 7.2921151467E-5;
            let dt_tx = dt_tx.to_seconds();
            let we = EARTH_OMEGA_E_WGS84 * dt_tx;
            let (we_cos, we_sin) = (we.cos(), we.sin());
            Matrix3::<f64>::new(
                we_cos, we_sin, 0.0_f64, -we_sin, we_cos, 0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64,
            )
        } else {
            Matrix3::<f64>::new(
                1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64, 0.0_f64, 0.0_f64, 0.0_f64, 1.0_f64,
            )
        };
        reworked.position = rot * interpolated.position;
        reworked
    }
    /*
     * Determine velocities
     */
    fn velocities(
        &self,
        t_tx: Epoch,
        sv: SV,
        interpolated: InterpolationResult,
    ) -> InterpolationResult {
        let mut reworked = interpolated;
        if let Some((p_ttx, p_pos)) = self.prev_sv_state.get(&sv) {
            let dt = (t_tx - *p_ttx).to_seconds();
            reworked.velocity = Some((interpolated.position - p_pos) / dt);
        }
        reworked
    }
    /*
     * Reworks solution a little bit
     */
    fn rework_solution(
        pvt: &mut PVTSolution,
        pvt_type: PVTSolutionType,
        apriori: &AprioriPosition,
        cfg: &Config,
    ) {
        if let Some(alt) = cfg.fixed_altitude {
            pvt.pos.z = apriori.ecef.z - alt;
            pvt.vel.z = 0.0_f64;
        }
        match pvt_type {
            PVTSolutionType::TimeOnly => {
                pvt.pos = Default::default();
                pvt.vel = Default::default();
            },
            _ => {},
        }
    }
}
