// SPP example (pseudo range based direct positioning).
// This is simply here to demonstrate how to operate the API, and does not generate actual results.
use gnss_rtk::prelude::{
    BaseStation as RTKBaseStation, Candidate, Carrier, ClockCorrection, Config, Duration, Epoch,
    Error, InvalidationCause, IonosphereBias, Method, Observation, OrbitalState,
    OrbitalStateProvider, Solver, TroposphereBias, SV,
};

// Orbit source example
struct Orbits {}

impl OrbitalStateProvider for Orbits {
    // For each requested "t" and "sv",
    // if we can, we should resolve the SV [OrbitalState].
    // If interpolation is to be used (depending on your apps), you can
    // use the interpolation order that we recommend here, or decide to ignore it.
    // If you're not in position to determine [OrbitalState], simply return None.
    // If None is returned for too long, this [Epoch] will eventually be dropped out,
    // and we will move on to the next
    fn next_at(&mut self, t: Epoch, sv: SV, order: usize) -> Option<OrbitalState> {
        let (x, y, z) = (0.0_f64, 0.0_f64, 0.0_f64);
        Some(OrbitalState::from_position((x, y, z)))
    }
}

// This example is direct positioning (not RTK), therefore
// the BaseStation returns Null all the time (== non existant)
struct BaseStation {}

impl RTKBaseStation for BaseStation {
    fn observe(&mut self, t: Epoch, sv: SV, carrier: Carrier) -> Option<Observation> {
        None // no differential positioning
    }
}

// Data source example
struct MyDataSource {}

impl MyDataSource {
    // Data source example
    fn new() -> Self {
        Self {}
    }
    // The objective here is to propose enough SV observations to resolve a solution.
    // Since our example only requires PseudoRange on a single frequency
    // we will limit ourselves to that.
    fn next(&mut self) -> Option<(Epoch, Vec<Candidate>)> {
        Some((
            // This must be the Epoch of observation,
            // ie, Sampling Instant
            Epoch::default(),
            vec![
                // Create a candidate from your Pseudo Range observation
                Candidate::new(
                    // Candidate Identity
                    SV::default(),
                    // Sampling Epoch. Must be identical for this grouping of candidates
                    Epoch::default(),
                    // You must provide a [ClockCorrection] for each candidate
                    ClockCorrection::without_relativistic_correction(Duration::from_nanoseconds(
                        100.0,
                    )),
                    // If you know the total group day for this Candidate, specify it here
                    None,
                    // List of observations
                    vec![Observation {
                        carrier: Carrier::L1, // example
                        pseudo: Some(3.0E6),  // example, this is raw observation
                        // Note that if you apply a min_snr preset,
                        // we might drop candidates that do not have this info
                        snr: None, // unknown
                        phase: None,
                        doppler: None,
                        ambiguity: None,
                    }],
                ),
                // Create all as many candidates as possible.
                // It's better to have more than needed, it leaves us more possibility in the election process.
            ],
        ))
    }
}

pub fn main() {
    // Build the Orbit source
    let orbits = Orbits {};

    // Build the Base station
    let base_station = BaseStation {};

    // The preset API is useful to quickly deploy depending on your application.
    // Static presets target static positioning.
    let cfg = Config::static_preset(Method::SPP); // Single Freq. Pseudo Range based

    // The API is pretty straightforward, it requires the Configuration preset to be
    // built ahead of time. The only difficulty is the design your data source and SV state provider. We interact with the SV state provider by means of a function pointer.
    let solver = Solver::new(
        &cfg,
        // We deploy without apriori knowledge.
        // The solver will initialize itself.
        None,
        // Tie the Orbit source
        orbits,
        // Tie the Base station
        Some(base_station),
    );

    // The solver needs to be mutable, due to the iteration process.
    let mut solver = solver.unwrap();

    let mut source = MyDataSource::new();

    // Browse your data source (This is an Example)
    while let Some((epoch, candidates)) = source.next() {
        // External bias sources are not very well supported yet.
        // We recommend you use the internal modeling.
        // This is done by specifying you simply have no knowledge
        // of the external biases:
        let ionod = IonosphereBias::default(); // Unknown
        let tropod = TroposphereBias::default(); // Unknown

        match solver.resolve(epoch, &candidates, &ionod, &tropod) {
            Ok((_epoch, solution)) => {
                // A solution was successfully resolved for this Epoch.
                // The position is expressed as absolute ECEF [m].
                let _position = solution.position;
                // The velocity vector is expressed as variations of absolute ECEF positions [m/s]
                let _velocity = solution.velocity;
                // Receiver offset to preset timescale
                let (_clock_offset, _timescale) = (solution.dt, solution.timescale);
                // More infos on SVs that contributed to this solution
                for info in solution.sv.values() {
                    // attitude
                    let (_el, _az) = (info.azimuth, info.elevation);
                    // Modeled (in this example) or simply copied ionosphere bias
                    // impacting selected signal from this spacecraft
                    let _bias_m = info.iono_bias;
                    // Modeled (in this example) or simply copied troposphere bias
                    // impacting selected signal from this spacecraft
                    let _bias_m = info.tropo_bias;
                    // Dilution of Precision informs on geometric performances
                    let (_tdop, _gdop, _pdop) = (solution.tdop, solution.gdop, solution.pdop);
                    // Determine the Vertical DoP for these lat,lon coordinates
                    let (lat_ddeg, lon_ddeg) = (45.0, 13.0);
                    let _vdop = solution.vdop(lat_ddeg, lon_ddeg);
                }
            },
            Err(Error::InvalidatedSolution(cause)) => match cause {
                InvalidationCause::FirstSolution => {
                    // The current behavior will always discard the first solution.
                    // We will always wind up here on the first iteration
                },
                _other => {
                    // Conditions are either degrading and no longer fit preset criteria,
                    // or we need more iterations to finally meet preset criteria.
                },
            },
            Err(_e) => {
                // Something went wrong, use "e" to get more info on what is wrong.
                // The most plausible cause is the lack of observations, at this point in time.
                // But that should not last for too long, otherwise something's wrong in your setup
            },
        }
    }
}
