use crate::prelude::{Candidate, Carrier, SPEED_OF_LIGHT_M_S};

/// Signal [Combination]
#[derive(Debug, Copy, Clone, Default)]
pub(crate) struct Combination {
    /// Lhs [Carrier] signal
    pub lhs: Carrier,

    /// Rhs reference [Carrier] signal
    pub rhs: Carrier,

    /// lambda
    pub lambda: f64,

    /// Value
    pub value: f64,
}

impl Candidate {
    /// Returns an IF code range [Combination] using prefered available frequencies
    pub(crate) fn code_if_combination(&self) -> Option<Combination> {
        let (c1, p1) = self.l1_pseudo_range()?;
        let f1_hz = c1.frequency_hz();

        let (cj, pj) = self.subsidary_pseudo_range()?;
        let fj_hz = cj.frequency_hz();

        let (f1pow, fjpow) = (f1_hz.powi(2), fj_hz.powi(2));

        let freq = f1_hz * fj_hz / (f1pow + fjpow).sqrt();

        Some(Combination {
            lhs: cj,
            rhs: c1,
            lambda: SPEED_OF_LIGHT_M_S / freq,
            value: (f1pow * p1 - fjpow * pj) / (f1pow - fjpow),
        })
    }

    /// Returns an IF phase range [Combination] using prefered available frequencies
    pub(crate) fn phase_if_combination(&self) -> Option<Combination> {
        let (c1, l1) = self.l1_phase_range()?;
        let f1_hz = c1.frequency_hz();

        let (cj, lj) = self.subsidary_phase_range()?;
        let fj_hz = cj.frequency_hz();

        let (f1pow, fjpow) = (f1_hz.powi(2), fj_hz.powi(2));

        let freq = f1_hz * fj_hz / (f1pow + fjpow).sqrt();

        Some(Combination {
            lhs: cj,
            rhs: c1,
            lambda: SPEED_OF_LIGHT_M_S / freq,
            value: (f1pow * l1 - fjpow * lj) / (f1pow - fjpow),
        })
    }

    /// Returns phase wide lane [Combination]
    pub(crate) fn phase_wl_combination(&self) -> Option<Combination> {
        let (c1, l1) = self.l1_phase_range()?;
        let f1 = c1.frequency_hz();

        let (cj, lj) = self.subsidary_phase_range()?;
        let fj = cj.frequency_hz();

        Some(Combination {
            lhs: cj,
            rhs: c1,
            lambda: SPEED_OF_LIGHT_M_S / (f1 - fj),
            value: (f1 * l1 - fj * lj) / (f1 - fj),
        })
    }

    /// Returns code narrow lane [Combination]
    pub(crate) fn code_nl_combination(&self) -> Option<Combination> {
        let (c1, p1) = self.l1_pseudo_range()?;
        let f1 = c1.frequency_hz();

        let (cj, pj) = self.subsidary_pseudo_range()?;
        let fj = cj.frequency_hz();

        Some(Combination {
            lhs: cj,
            rhs: c1,
            lambda: SPEED_OF_LIGHT_M_S / (f1 + fj),
            value: (f1 * p1 + fj * pj) / (f1 + fj),
        })
    }

    // /// Returns MW [Combination]
    // pub(crate) fn mw_combination(&self) -> Option<Combination> {
    //     let ph_w = self.phase_wl_combination()?;
    //     let pr_n = self.code_nl_combination()?;

    //     Some(Combination {
    //         lhs: ph_w.lhs,
    //         rhs: ph_w.rhs,
    //         value: ph_w.value - pr_n.value,
    //         lambda: ph_w.lambda,
    //     })
    // }

    /// Form GF [Combination]
    #[cfg(test)]
    pub(crate) fn code_gf_combination(&self) -> Option<Combination> {
        let (c1, p1) = self.l1_pseudo_range()?;
        let (c2, p2) = self.subsidary_pseudo_range()?;
        Some(Combination {
            lhs: c2,
            rhs: c1,
            value: p2 - p1,
            lambda: SPEED_OF_LIGHT_M_S / c1.frequency_hz(),
        })
    }

    /// Returns GF [Combination]
    pub(crate) fn phase_gf_combination(&self) -> Option<Combination> {
        let (c1, l1) = self.l1_phase_range()?;
        let (c2, l2) = self.subsidary_phase_range()?;
        Some(Combination {
            lhs: c2,
            rhs: c1,
            value: l1 - l2,
            lambda: SPEED_OF_LIGHT_M_S / c1.frequency_hz(),
        })
    }
}
