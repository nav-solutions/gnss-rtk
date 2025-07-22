/// [MeasurementSystem]s provide [Measurement]s that GNSS-RTK will
/// process to resolve P.V.T solutions.
pub struct Measurement {
    /// Satellite (signal source) as [SV]
    pub satellite: SV,

    /// [Carrier] is the sampled signal
    pub carrier: Carrier,

    /// Decoded pseudo-range, in meters.
    pub pseudo_range_m: Option<f64>,

    /// Phse range in meters.
    pub phase_range_m: Option<f64>,

    /// Doppler observation
    pub doppler: Option<f64>,
}

/// [MeasurementSystem] must be implemented by all measurement systems
/// forwarding measurements (oftentimes referred to as "observations").
/// When using differential navigation techniques, a minimum of two [MeasurementSystem]s are
/// required.
pub trait MeasurentSystem {
    /// Update your system and prepare for a new sampling [Epoch].
    /// GNSS-RTK follows the Timeserie you described
    fn new_epoch(&mut self, epoch: Epoch);

    /// Provide the solver with new measurements. 
    /// Return None when you cannot propose any more measurements for ongoing [Epoch].
    fn measure(&self) -> Option<Measurement>;
}
