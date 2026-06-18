#[derive(Debug, Copy, Clone)]
pub struct VectorContribution {
    /// Row #1 contribution
    pub row_1: f64,

    /// Row #2 contribution
    pub row_2: f64,

    /// Measurement standard deviation (m) for row #1.
    pub sigma_1: f64,

    /// Measurement standard deviation (m) for row #2.
    pub sigma_2: f64,
}

impl Default for VectorContribution {
    fn default() -> Self {
        // Unit deviations keep the weight matrix at identity (uniform
        // weighting) until a model fills them in.
        Self {
            row_1: 0.0,
            row_2: 0.0,
            sigma_1: 1.0,
            sigma_2: 1.0,
        }
    }
}
