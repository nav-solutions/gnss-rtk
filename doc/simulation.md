Simulation test bench
=====================

`src/tests/sim` contains a deterministic GNSS measurement simulator and a
scenario based test suite. It only exists under `cfg(test)`: the public API of
the library is not affected.

Why
---

Real data (RINEX) only tells us how far we are from a reference coordinate.
A simulator gives the *exact* truth of everything (position, receiver clock,
satellite clocks, atmosphere, integer ambiguities, cycle slips), so we can:

- prove that the processing chain is exact when the inputs are perfect
  (mm level, anything else is a bug),
- inject one error source at a time and check the solver reacts as physics
  predicts,
- measure performance metrics (RMS, 95th percentile, convergence time, clock
  accuracy, availability) on many scenarios, quickly, without network access.

The simulator already found several bugs (see the git history):
a 10x error on the gravitational constant, a wrong Shapiro formula,
the receiver clock leaking into the Sagnac correction, the ionosphere model
applied on top of the iono free combination, the postfit smoother feeding
back into the navigation filter, and the RTK far-field linearization.

Architecture
------------

| Module | Role |
|--------|------|
| `rng` | SplitMix64 + Box-Muller, bit exact reproducibility (no `rand` dependency) |
| `geodesy` | WGS84 helpers, independent from the solver (truth side) |
| `constellation` | Walker constellations (Galileo 24/3/1, GPS 24/6/1), ECEF positions |
| `scenario` | Pure description of a run: trajectory, clocks, noise, iono, tropo, slips, outages, base |
| `simulator` | Forward observation model, implements `OrbitSource`, `SpacebornBias`, `EnvironmentalBias`, `RTKBase` |
| `runner` | Feeds the public `Solver` epoch by epoch, compares to the truth, computes metrics |
| `suite` | The test suite (`model`, `absolute`, `rtk`, `kinematic`, `robustness`, `ppp`, `benchmark`) |

Forward model, for station r, satellite s, carrier f, receiver time tag T:

```text
t_rx = T - dt_r                                    (true reception time)
tau  = |R3(we.tau) x_s(t_rx - tau) - x_r(t_rx)| / c  (light time, Sagnac)
P_f  = c.tau + S + c.(dt_r - dt_s) + T_trop + I_f + e_P
L_f  = c.tau + S + c.(dt_r - dt_s) + T_trop - I_f + lambda_f.(N_f + slips) + e_L
```

- S: Shapiro delay, T_trop: Niel model (what the solver receives) + unmodeled ZWD,
  I_f = 40.3 STEC / f² (thin shell), the solver receives `(1 - model_error) * I`.
- Satellite clocks are quantized to `Duration` resolution so the truth and the
  solver inputs are identical.
- Noise: sigma(el) = sigma_zenith / sin(el).

Writing a scenario
------------------

```rust
let scenario = Scenario::realistic_static("my_scenario")   // geodetic noise, 20 TECu, TCXO
    .with_duration(3600.0, 30.0)
    .with_constellations(&[SimConstellation::GPS, SimConstellation::Galileo])
    .with_base_offset_enu(3_000.0, 4_000.0, 0.0)            // RTK base 5 km away
    .with_cycle_slip(CycleSlip { sv, at_s: 600.0, cycles_1: 1, cycles_2: 1 })
    .with_seed(42);

let sim = Simulator::new(scenario, earth_frame());

let report = run(&sim, &RunOptions::new(Method::CPP, Technique::Rtk)
    .with_params(UserParameters::new(UserProfile::Static, ClockProfile::Quartz)));

assert!(report.stats_from(600.0).rms_3d_m < 1.0);
assert!(report.convergence_time_s(1.0).unwrap() < 1800.0);
```

Running
-------

```bash
# whole suite (release is ~30x faster than debug for this crate)
cargo test --release

# only the simulation suite
cargo test --release sim::

# targets of the next milestones (expected to fail today)
cargo test --release sim:: -- --ignored

# performance report (markdown table), plus one CSV per run
GNSS_RTK_SIM_OUTPUT=target/sim cargo test --release performance_report -- --ignored --nocapture
```

Each CSV contains, per epoch: truth, estimate, ENU error, clock truth/estimate,
GDOP and the solver error if any.

Targets (`#[ignore]` tests)
---------------------------

- `kinematic::target_*`: a position-velocity state is needed for the
  pedestrian / car profiles to track moving receivers, and the default postfit
  smoother (q = 1E-6 m²) prevents any kinematic use.
- `ppp::target_*`: float PPP specifications (M2-M4).
