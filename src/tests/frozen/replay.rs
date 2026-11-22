//! Replays a frozen dataset through the [Solver] and builds a [RunReport].
use nalgebra::Vector3;

use crate::{
    prelude::{Config, Constellation, Method, Rc, Solver, UserParameters, UserProfile},
    tests::{
        TestEnvironment, almanac, earth_frame,
        ephemeris::NullEph,
        frozen::{Dataset, FrozenBase},
        navigation_test_config,
        sim::{
            geodesy::ecef_delta_to_enu,
            runner::{EpochOutcome, RunReport},
        },
        time::NullTime,
    },
};

/// Replay options
#[derive(Clone)]
pub struct Replay {
    pub dataset: Dataset,
    pub method: Method,
    pub rtk: bool,
    pub constellations: Vec<Constellation>,
    pub params: UserParameters,
    pub cfg: Config,
}

impl Replay {
    /// Static receiver, navigation filter evaluated without postfit smoother.
    pub fn new(
        dataset: Dataset,
        method: Method,
        rtk: bool,
        constellations: &[Constellation],
    ) -> Self {
        Self {
            dataset,
            method,
            rtk,
            constellations: constellations.to_vec(),
            params: UserParameters::new(UserProfile::Static, crate::prelude::ClockProfile::Quartz),
            cfg: navigation_test_config(method),
        }
    }

    pub fn with_params(mut self, params: UserParameters) -> Self {
        self.params = params;
        self
    }

    pub fn name(&self) -> String {
        let constellations = self
            .constellations
            .iter()
            .map(|c| match c {
                Constellation::GPS => "gps",
                Constellation::Galileo => "gal",
                Constellation::BeiDou => "bds",
                _ => "other",
            })
            .collect::<Vec<_>>()
            .join("+");

        format!(
            "frozen_{:?}_{}{}_{}",
            self.dataset,
            if self.rtk { "rtk-" } else { "" },
            self.method,
            constellations
        )
        .to_lowercase()
    }
}

/// Replays the dataset, returns the [RunReport].
///
/// Truth: rover antenna reference point, in ITRF2014 for the absolute techniques
/// (frame of the precise products), in ETRF2000 (frame of the base coordinates) for RTK.
pub fn replay(opts: &Replay) -> RunReport {
    let data = opts.dataset.load();
    let frame = earth_frame();

    let mut solver = Solver::new_survey(
        almanac(),
        frame,
        opts.cfg.clone(),
        Rc::new(NullEph {}),
        Rc::new(&data.orbits),
        Rc::new(&data.clocks),
        Rc::new(TestEnvironment::new()),
        NullTime {},
    );

    let base = FrozenBase {
        station: &data.base,
        constellations: opts.constellations.clone(),
    };

    let truth = if opts.rtk {
        data.rover.arp_ecef_m()
    } else {
        data.rover.itrf2014_arp_ecef_m()
    };
    let mut outcomes = Vec::with_capacity(data.rover.epochs.len());

    for (k, epoch) in data.rover.epochs.keys().enumerate() {
        let candidates = data.rover.candidates(*epoch, &opts.constellations);

        let result = if opts.rtk {
            solver.rtk(*epoch, opts.params, &candidates, &base)
        } else {
            solver.ppp(*epoch, opts.params, &candidates)
        };

        let mut outcome = EpochOutcome {
            k,
            t_s: (*epoch - data.rover.t0).to_seconds(),
            num_candidates: candidates.len(),
            truth_m: truth,
            truth_clock_s: f64::NAN,
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
            Err(e) => outcome.failure = Some(e.to_string()),
        }

        outcomes.push(outcome);
    }

    let report = RunReport {
        name: opts.name(),
        outcomes,
        orbit_requests: 0,
    };

    log::info!("{report}");
    report.dump_if_requested();
    report
}
