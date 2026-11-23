//! Deterministic random number generation for the simulator.
//!
//! We do not depend on `rand` on purpose: scenarios must be bit-exact
//! reproducible across `rand` releases and platforms.

/// SplitMix64 generator with Gaussian sampling (Box-Muller).
#[derive(Debug, Clone)]
pub struct SimRng {
    state: u64,
    spare: Option<f64>,
}

impl SimRng {
    pub fn new(seed: u64) -> Self {
        Self {
            state: seed ^ 0x9E37_79B9_7F4A_7C15,
            spare: None,
        }
    }

    /// Derive an independent stream (for one satellite, one signal...).
    pub fn fork(&self, salt: u64) -> Self {
        let mut s = Self::new(
            self.state
                .wrapping_add(salt.wrapping_mul(0xD1B5_4A32_D192_ED03)),
        );
        // decorrelate
        for _ in 0..4 {
            s.next_u64();
        }
        s
    }

    pub fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    /// Uniform in [0, 1)
    pub fn uniform(&mut self) -> f64 {
        ((self.next_u64() >> 11) as f64) / ((1u64 << 53) as f64)
    }

    /// Uniform in [lo, hi)
    pub fn uniform_in(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.uniform()
    }

    /// Standard normal sample
    pub fn gaussian(&mut self) -> f64 {
        if let Some(v) = self.spare.take() {
            return v;
        }

        let u1 = self.uniform().max(1.0E-300);
        let u2 = self.uniform();
        let r = (-2.0 * u1.ln()).sqrt();
        let theta = 2.0 * std::f64::consts::PI * u2;

        self.spare = Some(r * theta.sin());
        r * theta.cos()
    }

    /// Normal sample N(0, sigma²)
    pub fn normal(&mut self, sigma: f64) -> f64 {
        sigma * self.gaussian()
    }
}

#[cfg(test)]
mod test {
    use super::SimRng;

    #[test]
    fn reproducible() {
        let mut a = SimRng::new(1);
        let mut b = SimRng::new(1);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn gaussian_moments() {
        let mut rng = SimRng::new(123);
        let n = 200_000;
        let samples: Vec<f64> = (0..n).map(|_| rng.gaussian()).collect();
        let mean = samples.iter().sum::<f64>() / n as f64;
        let var = samples.iter().map(|v| (v - mean).powi(2)).sum::<f64>() / n as f64;
        assert!(mean.abs() < 0.01, "mean={mean}");
        assert!((var - 1.0).abs() < 0.02, "var={var}");
    }
}
