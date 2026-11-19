//! Scenario runner and performance metrics.
//!
//! A run feeds every simulated epoch to a [Solver] through its public API
//! and compares each [PVTSolution] to the truth. The resulting [RunReport]
//! exposes the classical GNSS performance metrics.
//!
//! Set `GNSS_RTK_SIM_OUTPUT=<directory>` to dump every run as CSV
//! (one file per scenario) for plotting / regression tracking.
use std::fmt;

use nalgebra::Vector3;

use crate::{
    prelude::{Config, Method, Rc, Solver, UserParameters},
    tests::{
        almanac,
        ephemeris::NullEph,
        sim::{geodesy::ecef_delta_to_enu, simulator::Simulator},
        time::NullTime,
    },
};

/// Navigation technique
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Technique {
    /// Absolute navigation (`Solver::ppp`)
    Absolute,
    /// Differential navigation (`Solver::rtk`)
    Rtk,
}

/// How the solver is initialized
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Initialization {
    /// No apriori: survey from scratch (Bancroft)
    Survey,
    /// Apriori given with the following ECEF offset from the truth (m)
    Apriori(f64, f64, f64),
}

/// Run options
#[derive(Debug, Clone)]
pub struct RunOptions {
    pub cfg: Config,
    pub params: UserParameters,
    pub technique: Technique,
    pub initialization: Initialization,
}

impl RunOptions {
    pub fn new(method: Method, technique: Technique) -> Self {
        let mut cfg = Config::default().with_navigation_method(method);
        // the synthetic orbits are expressed in the earth fixed frame directly:
        // the eclipse model would randomly drop satellites.
        cfg.max_eclipse_rate_percent = 100.0;

        Self {
            cfg,
            params: UserParameters::default(),
            technique,
            initialization: Initialization::Survey,
        }
    }

    pub fn with_params(mut self, params: UserParameters) -> Self {
        self.params = params;
        self
    }

    pub fn with_initialization(mut self, init: Initialization) -> Self {
        self.initialization = init;
        self
    }

    pub fn with_config(mut self, f: impl FnOnce(&mut Config)) -> Self {
        f(&mut self.cfg);
        self
    }
}

/// Outcome of one epoch
#[derive(Debug, Clone)]
pub struct EpochOutcome {
    pub k: usize,
    pub t_s: f64,
    pub num_candidates: usize,
    pub truth_m: Vector3<f64>,
    pub truth_clock_s: f64,
    /// Some(estimated position, estimated clock) when the solver produced a solution
    pub solution: Option<(Vector3<f64>, f64)>,
    /// ENU error (m)
    pub error_enu_m: Option<Vector3<f64>>,
    pub gdop: Option<f64>,
    pub tdop: Option<f64>,
    pub failure: Option<String>,
}

impl EpochOutcome {
    pub fn error_3d_m(&self) -> Option<f64> {
        self.error_enu_m.map(|e| e.norm())
    }

    pub fn error_h_m(&self) -> Option<f64> {
        self.error_enu_m.map(|e| (e[0] * e[0] + e[1] * e[1]).sqrt())
    }

    pub fn error_v_m(&self) -> Option<f64> {
        self.error_enu_m.map(|e| e[2].abs())
    }

    /// Clock error (s)
    pub fn clock_error_s(&self) -> Option<f64> {
        self.solution.map(|(_, clk)| clk - self.truth_clock_s)
    }
}

/// Summary statistics over a set of epochs
#[derive(Debug, Clone, Copy, Default)]
pub struct Stats {
    pub count: usize,
    pub rms_e_m: f64,
    pub rms_n_m: f64,
    pub rms_u_m: f64,
    pub rms_h_m: f64,
    pub rms_3d_m: f64,
    pub p95_3d_m: f64,
    pub max_3d_m: f64,
    pub mean_enu_m: (f64, f64, f64),
}

/// Worst values (per epoch) over a set of epochs
#[derive(Debug, Clone, Copy, Default)]
pub struct Worst {
    /// Number of solutions considered
    pub count: usize,
    /// max |East|, |North|, |Up| errors (m)
    pub e_m: f64,
    pub n_m: f64,
    pub u_m: f64,
    /// max |clock error| (s), when the truth is known
    pub dt_s: Option<f64>,
    /// max deviation of the clock estimate from its linear fit (s):
    /// clock estimation noise, usable without truth (steered receivers).
    pub dt_jitter_s: f64,
    pub gdop: f64,
    pub tdop: f64,
}

/// Result of a full run
#[derive(Debug, Clone)]
pub struct RunReport {
    pub name: String,
    pub outcomes: Vec<EpochOutcome>,
    pub orbit_requests: usize,
}

fn rms(values: impl Iterator<Item = f64>) -> (f64, usize) {
    let mut sum = 0.0;
    let mut n = 0;
    for v in values {
        sum += v * v;
        n += 1;
    }
    if n == 0 {
        (f64::NAN, 0)
    } else {
        ((sum / n as f64).sqrt(), n)
    }
}

fn percentile(mut values: Vec<f64>, p: f64) -> f64 {
    if values.is_empty() {
        return f64::NAN;
    }
    values.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let idx = ((p / 100.0) * (values.len() - 1) as f64).round() as usize;
    values[idx.min(values.len() - 1)]
}

impl RunReport {
    pub fn num_epochs(&self) -> usize {
        self.outcomes.len()
    }

    pub fn num_solutions(&self) -> usize {
        self.outcomes
            .iter()
            .filter(|o| o.solution.is_some())
            .count()
    }

    /// Ratio of epochs that produced a solution
    pub fn availability(&self) -> f64 {
        if self.outcomes.is_empty() {
            return 0.0;
        }
        self.num_solutions() as f64 / self.num_epochs() as f64
    }

    /// Statistics over solutions obtained at `t_s >= from_s`
    pub fn stats_from(&self, from_s: f64) -> Stats {
        let errors: Vec<Vector3<f64>> = self
            .outcomes
            .iter()
            .filter(|o| o.t_s >= from_s)
            .filter_map(|o| o.error_enu_m)
            .collect();

        if errors.is_empty() {
            return Stats::default();
        }

        let n = errors.len() as f64;
        let (rms_e_m, count) = rms(errors.iter().map(|e| e[0]));
        let (rms_n_m, _) = rms(errors.iter().map(|e| e[1]));
        let (rms_u_m, _) = rms(errors.iter().map(|e| e[2]));
        let (rms_h_m, _) = rms(errors.iter().map(|e| (e[0] * e[0] + e[1] * e[1]).sqrt()));
        let (rms_3d_m, _) = rms(errors.iter().map(|e| e.norm()));
        let norms: Vec<f64> = errors.iter().map(|e| e.norm()).collect();

        Stats {
            count,
            rms_e_m,
            rms_n_m,
            rms_u_m,
            rms_h_m,
            rms_3d_m,
            p95_3d_m: percentile(norms.clone(), 95.0),
            max_3d_m: norms.iter().cloned().fold(0.0, f64::max),
            mean_enu_m: (
                errors.iter().map(|e| e[0]).sum::<f64>() / n,
                errors.iter().map(|e| e[1]).sum::<f64>() / n,
                errors.iter().map(|e| e[2]).sum::<f64>() / n,
            ),
        }
    }

    /// Worst values over solutions obtained at `t_s >= from_s`
    pub fn worst_from(&self, from_s: f64) -> Worst {
        let sols = self
            .outcomes
            .iter()
            .filter(|o| o.t_s >= from_s && o.solution.is_some())
            .collect::<Vec<_>>();

        let mut w = Worst {
            count: sols.len(),
            ..Default::default()
        };

        if sols.is_empty() {
            return w;
        }

        for o in sols.iter() {
            if let Some(e) = o.error_enu_m {
                w.e_m = w.e_m.max(e[0].abs());
                w.n_m = w.n_m.max(e[1].abs());
                w.u_m = w.u_m.max(e[2].abs());
            }
            if let Some(dt) = o.clock_error_s() {
                if dt.is_finite() {
                    w.dt_s = Some(w.dt_s.unwrap_or(0.0).max(dt.abs()));
                }
            }
            w.gdop = w.gdop.max(o.gdop.unwrap_or(0.0));
            w.tdop = w.tdop.max(o.tdop.unwrap_or(0.0));
        }

        // clock estimation noise: deviation from a linear fit
        let pts = sols
            .iter()
            .map(|o| (o.t_s, o.solution.unwrap().1))
            .collect::<Vec<_>>();

        let n = pts.len() as f64;
        let mean_t = pts.iter().map(|p| p.0).sum::<f64>() / n;
        let mean_c = pts.iter().map(|p| p.1).sum::<f64>() / n;
        let var_t = pts.iter().map(|p| (p.0 - mean_t).powi(2)).sum::<f64>();
        let slope = if var_t > 0.0 {
            pts.iter()
                .map(|p| (p.0 - mean_t) * (p.1 - mean_c))
                .sum::<f64>()
                / var_t
        } else {
            0.0
        };

        w.dt_jitter_s = pts
            .iter()
            .map(|p| (p.1 - (mean_c + slope * (p.0 - mean_t))).abs())
            .fold(0.0, f64::max);

        w
    }

    /// Statistics over the whole run
    pub fn stats(&self) -> Stats {
        self.stats_from(f64::NEG_INFINITY)
    }

    /// Convergence time: first instant after which the 3D error stays
    /// below `threshold_m` until the end of the run. None if never converged.
    pub fn convergence_time_s(&self, threshold_m: f64) -> Option<f64> {
        let mut candidate = None;
        for o in self.outcomes.iter() {
            match o.error_3d_m() {
                Some(err) if err < threshold_m => {
                    if candidate.is_none() {
                        candidate = Some(o.t_s);
                    }
                },
                _ => candidate = None,
            }
        }
        candidate
    }

    /// RMS of the receiver clock error (ns), from `from_s`
    pub fn clock_rms_ns_from(&self, from_s: f64) -> f64 {
        rms(self
            .outcomes
            .iter()
            .filter(|o| o.t_s >= from_s)
            .filter_map(|o| o.clock_error_s())
            .map(|e| e * 1.0E9))
        .0
    }

    /// First solver failures (for diagnostics)
    pub fn failures(&self) -> Vec<(usize, String)> {
        self.outcomes
            .iter()
            .filter_map(|o| o.failure.clone().map(|f| (o.k, f)))
            .collect()
    }

    /// CSV dump
    pub fn to_csv(&self) -> String {
        let mut s = String::from(
            "k,t_s,nsat,truth_x_m,truth_y_m,truth_z_m,est_x_m,est_y_m,est_z_m,err_e_m,err_n_m,err_u_m,err_3d_m,truth_clk_s,est_clk_s,gdop,failure\n",
        );
        for o in self.outcomes.iter() {
            let (est, clk) = match o.solution {
                Some((p, c)) => (
                    format!("{:.4},{:.4},{:.4}", p[0], p[1], p[2]),
                    format!("{:.12e}", c),
                ),
                None => (",,".to_string(), String::new()),
            };
            let err = match o.error_enu_m {
                Some(e) => format!("{:.4},{:.4},{:.4},{:.4}", e[0], e[1], e[2], e.norm()),
                None => ",,,".to_string(),
            };
            s.push_str(&format!(
                "{},{},{},{:.4},{:.4},{:.4},{},{},{:.12e},{},{},{}\n",
                o.k,
                o.t_s,
                o.num_candidates,
                o.truth_m[0],
                o.truth_m[1],
                o.truth_m[2],
                est,
                err,
                o.truth_clock_s,
                clk,
                o.gdop.map(|g| format!("{g:.3}")).unwrap_or_default(),
                o.failure.clone().unwrap_or_default().replace(',', ";"),
            ));
        }
        s
    }

    /// Writes the CSV dump when `GNSS_RTK_SIM_OUTPUT` is defined.
    pub fn dump_if_requested(&self) {
        if let Ok(dir) = std::env::var("GNSS_RTK_SIM_OUTPUT") {
            let _ = std::fs::create_dir_all(&dir);
            let path = std::path::Path::new(&dir).join(format!("{}.csv", self.name));
            if let Err(e) = std::fs::write(&path, self.to_csv()) {
                log::error!("failed to write {}: {e}", path.display());
            }
        }
    }
}

impl fmt::Display for RunReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = self.stats();
        write!(
            f,
            "{:<40} avail={:>5.1}% rms(E,N,U)=({:.3},{:.3},{:.3})m rms3d={:.3}m p95={:.3}m max={:.3}m clk_rms={:.2}ns",
            self.name,
            self.availability() * 100.0,
            s.rms_e_m,
            s.rms_n_m,
            s.rms_u_m,
            s.rms_3d_m,
            s.p95_3d_m,
            s.max_3d_m,
            self.clock_rms_ns_from(f64::NEG_INFINITY),
        )
    }
}

/// Runs a simulated scenario through the [Solver].
pub fn run(sim: &Simulator, opts: &RunOptions) -> RunReport {
    let scenario = sim.scenario();

    let apriori = match opts.initialization {
        Initialization::Survey => None,
        Initialization::Apriori(dx, dy, dz) => {
            let r0 = scenario.rover.initial_position_m();
            Some((r0[0] + dx, r0[1] + dy, r0[2] + dz))
        },
    };

    let mut solver = Solver::new(
        almanac(),
        sim.frame(),
        opts.cfg.clone(),
        Rc::new(NullEph {}),
        Rc::new(sim.orbit_source()),
        Rc::new(sim.spaceborn_biases()),
        Rc::new(sim.environment()),
        NullTime {},
        apriori,
    );

    let base = if opts.technique == Technique::Rtk {
        Some(sim.rtk_base())
    } else {
        None
    };

    let mut outcomes = Vec::with_capacity(scenario.num_epochs());

    for sim_epoch in sim.rover_epochs() {
        let result = match &base {
            Some(base) => solver.rtk(sim_epoch.epoch, opts.params, &sim_epoch.candidates, base),
            None => solver.ppp(sim_epoch.epoch, opts.params, &sim_epoch.candidates),
        };

        let truth = sim_epoch.truth.position_m;

        let mut outcome = EpochOutcome {
            k: sim_epoch.k,
            t_s: sim_epoch.t_s,
            num_candidates: sim_epoch.candidates.len(),
            truth_m: truth,
            truth_clock_s: sim_epoch.truth.clock_offset_s,
            solution: None,
            error_enu_m: None,
            gdop: None,
            tdop: None,
            failure: None,
        };

        match result {
            Ok(pvt) => {
                let est = Vector3::new(pvt.pos_m.0, pvt.pos_m.1, pvt.pos_m.2);
                outcome.solution = Some((est, pvt.clock_offset_s));
                outcome.error_enu_m = Some(ecef_delta_to_enu(&truth, &(est - truth)));
                outcome.gdop = Some(pvt.gdop);
                outcome.tdop = Some(pvt.tdop);
            },
            Err(e) => {
                outcome.failure = Some(e.to_string());
            },
        }

        outcomes.push(outcome);
    }

    let report = RunReport {
        name: scenario.name.clone(),
        outcomes,
        orbit_requests: sim.orbit_requests(),
    };

    log::info!("{report}");
    report.dump_if_requested();
    report
}
