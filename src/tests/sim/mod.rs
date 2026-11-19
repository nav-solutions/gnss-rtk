//! GNSS measurement simulator and scenario based test bench.
//!
//! Layout:
//! - [rng]: deterministic random generation (reproducible runs)
//! - [geodesy]: independent WGS84 helpers (truth side)
//! - [constellation]: synthetic Walker constellations
//! - [scenario]: pure scenario descriptions + presets
//! - [simulator]: forward observation model, implements all solver input traits
//! - [runner]: runs a scenario through the public `Solver` API and computes metrics
//! - [suite]: the scenario based test suite
//!
//! Everything lives under `cfg(test)`: the public API of the library is untouched.
//! See `doc/simulation.md`.

// toolkit API: not every helper is used by the current suite
#![allow(dead_code)]
pub mod constellation;
pub mod geodesy;
pub mod rng;
pub mod runner;
pub mod scenario;
pub mod simulator;

mod suite;

pub use constellation::SimConstellation;
pub use runner::{Initialization, RunOptions, RunReport, Technique, run};
pub use scenario::{
    ClockModel, CycleSlip, IonoModel, NoiseModel, Outage, Scenario, Trajectory, TropoModel,
};
pub use simulator::{Simulator, Station};
