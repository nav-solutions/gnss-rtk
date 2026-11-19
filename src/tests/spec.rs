//! Performance contract of each solver preset.
//!
//! Every test that runs the navigation solver in nominal conditions checks
//! its solutions against the contract of its preset (also published in the README):
//!
//! | Preset   | max abs(East), abs(North) | max abs(Up) | clock (truth / stability) | max GDOP | max TDOP |
//! |----------|----------|--------|---------------|-----|------------------|
//! | SPP      | 5 m      | 10 m   | 50 ns / 30 ns | 5   | 3                |
//! | CPP      | 1 m      | 2 m    | 10 ns / 10 ns | 5   | 3                |
//! | PPP      | 0.3 m    | 0.5 m  | 2 ns / 2 ns   | 5   | 3                |
//! | RTK-SPP  | 0.5 m    | 0.5 m  | n/a           | 5   | 0 (not estimated)|
//! | RTK-CPP  | 0.5 m    | 0.5 m  | n/a           | 5   | 0 (not estimated)|
//! | RTK-PPP  | 0.5 m    | 0.5 m  | n/a           | 5   | 0 (not estimated)|
//!
//! - Errors are expressed in the local East / North / Up frame of the truth:
//!   ECEF axes mix horizontal and vertical errors depending on the latitude.
//! - Every epoch, after the convergence period of the scenario, must comply.
//! - Clock: against the truth when known (simulation), otherwise the stability of
//!   the estimate (max deviation from its linear fit, steered receivers).
//! - Nominal conditions: static or low dynamics receiver, open sky, GDOP <= 5,
//!   ionosphere model for SPP, RTK baseline up to 65 km.
use nalgebra::Vector3;

use crate::{
    prelude::{Method, PVTSolution},
    tests::sim::{geodesy::ecef_delta_to_enu, runner::RunReport},
};

/// Solver preset
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Preset {
    SPP,
    CPP,
    PPP,
    RtkSPP,
    RtkCPP,
    RtkPPP,
}

/// Performance contract of a [Preset]
#[derive(Debug, Clone, Copy)]
pub struct Contract {
    /// max |East| and |North| error (m)
    pub max_horizontal_axis_m: f64,
    /// max |Up| error (m)
    pub max_vertical_m: f64,
    /// max |clock error| (s), when the truth is known (absolute navigation only)
    pub max_clock_error_s: Option<f64>,
    /// max deviation of the clock estimate from its linear fit (s)
    pub max_clock_jitter_s: Option<f64>,
    pub max_gdop: f64,
    /// None: the clock is not estimated (TDOP must be null)
    pub max_tdop: Option<f64>,
}

impl Preset {
    pub fn new(method: Method, rtk: bool) -> Self {
        match (method, rtk) {
            (Method::SPP, false) => Self::SPP,
            (Method::CPP, false) => Self::CPP,
            (Method::PPP, false) => Self::PPP,
            (Method::SPP, true) => Self::RtkSPP,
            (Method::CPP, true) => Self::RtkCPP,
            (Method::PPP, true) => Self::RtkPPP,
        }
    }

    pub fn contract(&self) -> Contract {
        match self {
            Self::SPP => Contract {
                max_horizontal_axis_m: 5.0,
                max_vertical_m: 10.0,
                max_clock_error_s: Some(50.0E-9),
                max_clock_jitter_s: Some(30.0E-9),
                max_gdop: 5.0,
                max_tdop: Some(3.0),
            },
            Self::CPP => Contract {
                max_horizontal_axis_m: 1.0,
                max_vertical_m: 2.0,
                max_clock_error_s: Some(10.0E-9),
                max_clock_jitter_s: Some(10.0E-9),
                max_gdop: 5.0,
                max_tdop: Some(3.0),
            },
            Self::PPP => Contract {
                max_horizontal_axis_m: 0.3,
                max_vertical_m: 0.5,
                max_clock_error_s: Some(2.0E-9),
                max_clock_jitter_s: Some(2.0E-9),
                max_gdop: 5.0,
                max_tdop: Some(3.0),
            },
            Self::RtkSPP | Self::RtkCPP | Self::RtkPPP => Contract {
                max_horizontal_axis_m: 0.5,
                max_vertical_m: 0.5,
                max_clock_error_s: None,
                max_clock_jitter_s: None,
                max_gdop: 5.0,
                max_tdop: None,
            },
        }
    }
}

impl Contract {
    /// Violations of a single solution (empty when compliant)
    fn epoch_violations(
        &self,
        enu_m: &Vector3<f64>,
        clock_error_s: Option<f64>,
        gdop: f64,
        tdop: f64,
    ) -> Vec<String> {
        let mut v = Vec::new();

        for (axis, err) in [("East", enu_m[0]), ("North", enu_m[1])] {
            if !(err.abs() <= self.max_horizontal_axis_m) {
                v.push(format!(
                    "|{axis}| = {:.3} m > {:.3} m",
                    err.abs(),
                    self.max_horizontal_axis_m
                ));
            }
        }

        if !(enu_m[2].abs() <= self.max_vertical_m) {
            v.push(format!(
                "|Up| = {:.3} m > {:.3} m",
                enu_m[2].abs(),
                self.max_vertical_m
            ));
        }

        if let (Some(max), Some(err)) = (self.max_clock_error_s, clock_error_s) {
            if err.is_finite() && err.abs() > max {
                v.push(format!(
                    "|dt| = {:.2} ns > {:.2} ns",
                    err.abs() * 1E9,
                    max * 1E9
                ));
            }
        }

        if !(gdop <= self.max_gdop) {
            v.push(format!("GDOP = {gdop:.2} > {:.2}", self.max_gdop));
        }

        match self.max_tdop {
            Some(max) => {
                if !(tdop <= max) {
                    v.push(format!("TDOP = {tdop:.2} > {max:.2}"));
                }
            },
            None => {
                if tdop != 0.0 {
                    v.push(format!(
                        "TDOP = {tdop:.2}, expected 0 (clock not estimated)"
                    ));
                }
            },
        }

        v
    }

    /// Violations of a run, for solutions obtained at `t_s >= from_s`
    pub fn report_violations(&self, report: &RunReport, from_s: f64) -> Vec<String> {
        let mut v = Vec::new();

        for o in report.outcomes.iter().filter(|o| o.t_s >= from_s) {
            if let (Some(enu), Some(gdop), Some(tdop)) = (o.error_enu_m, o.gdop, o.tdop) {
                for violation in self.epoch_violations(&enu, o.clock_error_s(), gdop, tdop) {
                    v.push(format!("t={}s: {violation}", o.t_s));
                }
            }
        }

        // clock stability, only when the clock truth is unknown (real data):
        // a free running clock legitimately departs from a linear model.
        if let Some(max) = self.max_clock_jitter_s {
            let worst = report.worst_from(from_s);
            if worst.dt_s.is_none() && worst.count > 2 && worst.dt_jitter_s > max {
                v.push(format!(
                    "clock jitter = {:.2} ns > {:.2} ns",
                    worst.dt_jitter_s * 1E9,
                    max * 1E9
                ));
            }
        }

        v
    }
}

/// Asserts that every solution of the run (after `from_s`) complies with the
/// contract of `preset`. Lists the first violations otherwise.
pub fn assert_contract(preset: Preset, report: &RunReport, from_s: f64) {
    let violations = preset.contract().report_violations(report, from_s);

    assert!(
        report.num_solutions() > 0,
        "{}: no solution to check against the {preset:?} contract",
        report.name
    );

    assert!(
        violations.is_empty(),
        "{}: {} violation(s) of the {preset:?} contract, first ones:\n{}",
        report.name,
        violations.len(),
        violations
            .iter()
            .take(8)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// Asserts that a single [PVTSolution] complies with the contract of `preset`
/// (truth: ECEF coordinates in meters, clock truth unknown).
pub fn assert_pvt_contract(preset: Preset, pvt: &PVTSolution, truth_ecef_m: (f64, f64, f64)) {
    let truth = Vector3::new(truth_ecef_m.0, truth_ecef_m.1, truth_ecef_m.2);
    let est = Vector3::new(pvt.pos_m.0, pvt.pos_m.1, pvt.pos_m.2);
    let enu = ecef_delta_to_enu(&truth, &(est - truth));

    let violations = preset
        .contract()
        .epoch_violations(&enu, None, pvt.gdop, pvt.tdop);

    assert!(
        violations.is_empty(),
        "{}: violation(s) of the {preset:?} contract:\n{}",
        pvt.epoch,
        violations.join("\n")
    );
}

/// Same as [assert_pvt_contract], for a case with a known derogation:
/// violations are reported (warning) but do not fail the test.
pub fn check_pvt_contract_derogated(
    preset: Preset,
    pvt: &PVTSolution,
    truth_ecef_m: (f64, f64, f64),
    derogation: &str,
) {
    let truth = Vector3::new(truth_ecef_m.0, truth_ecef_m.1, truth_ecef_m.2);
    let est = Vector3::new(pvt.pos_m.0, pvt.pos_m.1, pvt.pos_m.2);
    let enu = ecef_delta_to_enu(&truth, &(est - truth));

    for violation in preset
        .contract()
        .epoch_violations(&enu, None, pvt.gdop, pvt.tdop)
    {
        log::warn!(
            "{}: {preset:?} contract derogation ({derogation}): {violation}",
            pvt.epoch
        );
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn contract_violations() {
        let c = Preset::CPP.contract();
        assert!(
            c.epoch_violations(&Vector3::new(0.5, -0.9, 1.9), Some(5E-9), 2.0, 1.0)
                .is_empty()
        );
        assert_eq!(
            c.epoch_violations(&Vector3::new(1.1, 0.0, 0.0), None, 2.0, 1.0)
                .len(),
            1
        );
        assert_eq!(
            c.epoch_violations(&Vector3::new(0.0, 0.0, -2.5), None, 2.0, 1.0)
                .len(),
            1
        );
        assert_eq!(
            c.epoch_violations(&Vector3::zeros(), Some(-20E-9), 6.0, 4.0)
                .len(),
            3
        );
        assert_eq!(
            c.epoch_violations(&Vector3::new(f64::NAN, 0.0, 0.0), None, 1.0, 1.0)
                .len(),
            1
        );

        let rtk = Preset::RtkCPP.contract();
        assert!(
            rtk.epoch_violations(&Vector3::new(0.4, 0.4, 0.4), None, 2.0, 0.0)
                .is_empty()
        );
        assert_eq!(
            rtk.epoch_violations(&Vector3::zeros(), None, 2.0, 0.5)
                .len(),
            1
        );
    }

    #[test]
    fn contracts_are_ordered() {
        let spp = Preset::SPP.contract();
        let cpp = Preset::CPP.contract();
        let ppp = Preset::PPP.contract();
        assert!(spp.max_horizontal_axis_m > cpp.max_horizontal_axis_m);
        assert!(cpp.max_horizontal_axis_m > ppp.max_horizontal_axis_m);
        assert!(spp.max_vertical_m > cpp.max_vertical_m);
        assert!(cpp.max_vertical_m > ppp.max_vertical_m);
    }
}
