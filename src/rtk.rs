use crate::prelude::{Candidate, Epoch};

/// Any [RTKBase] provides remote data by implementing the [RemoteSource] trait.
pub trait RTKBase {
    /// Provide a meaningful name of each reference station.
    /// This is only used when reporting status and results and is not vital.
    fn name(&self) -> String;

    /// Notifies the remote station that we're about to process a new [Epoch].
    /// Since GNSS-RTK expects measurements in chronological order, and simply
    /// follows the measurements, you should consider that any past [Epoch]
    /// will no longer be processed.
    fn new_epoch(&mut self, epoch: Epoch);

    /// Provide as many remote synchronous measurements (wrapped as [Candidate]s), at specified sampling [Epoch], as you can.
    /// In practice, it is most often impossible to be truly be synchronous, you should
    /// always minimize the time-difference, especially when the base is moving.
    fn observe(&self, epoch: Epoch) -> Vec<Candidate>;

    /// Provide the reference position of the base station, at sampling [Epoch].
    /// Any moving base station should keep its position up to date and garatee the age of the
    /// observations are kept to a minimal, the time-difference between directly related to the
    /// base velocity.
    /// This information is mandatory. Any error on this value will automatically
    /// decrease the accuracy of the solution.
    fn reference_position_ecef_m(&self, epoch: Epoch) -> (f64, f64, f64);
}

pub(crate) struct NullRTK {}

impl RTKBase for NullRTK {
    fn name(&self) -> String {
        "NULL".to_string()
    }

    fn new_epoch(&mut self, _: Epoch) {}

    fn observe(&self, _: Epoch) -> Vec<Candidate> {
        Default::default()
    }

    fn reference_position_ecef_m(&self, _: Epoch) -> (f64, f64, f64) {
        Default::default()
    }
}
