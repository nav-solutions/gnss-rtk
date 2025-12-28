use log::debug;

use std::collections::HashMap;

use crate::{
    constants::SPEED_OF_LIGHT_M_S,
    prelude::{Carrier, SV},
};

#[derive(Default, Debug)]
pub(crate) struct Difference {
    /// Code difference
    pub code: Option<(Carrier, f64)>,

    /// Extra code difference
    pub code_j: Option<(Carrier, f64)>,

    /// Code IF difference
    pub code_if: Option<(Carrier, f64)>,

    /// Phase difference
    pub phase: Option<(Carrier, f64)>,

    /// Extra phase difference
    pub phase_j: Option<(Carrier, f64)>,

    /// Phase IF difference
    pub phase_if: Option<(Carrier, f64, f64)>,

    /// Lw difference
    pub lw: Option<(Carrier, f64, f64)>,

    /// Cn difference
    pub cn: Option<(Carrier, f64, f64)>,
}

impl std::fmt::Display for Difference {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.is_some() {
            writeln!(f)?;
        }

        if let Some((carrier, code)) = &self.code {
            writeln!(f, "code({carrier})={code}")?;
        }

        if let Some((carrier, code)) = &self.code_j {
            writeln!(f, "code({carrier})={code}")?;
        }

        if let Some((carrier, code)) = &self.code_if {
            writeln!(f, "code_if({carrier})={code}")?;
        }

        if let Some((carrier, phase)) = &self.phase {
            writeln!(f, "phase({carrier})={phase}")?;
        }

        if let Some((carrier, phase)) = &self.phase_j {
            writeln!(f, "phase({carrier})={phase}")?;
        }

        if let Some((carrier, _, phase_if)) = &self.phase_if {
            writeln!(f, "phase_if({carrier})={phase_if} ")?;
        }

        if let Some((carrier, _, lw)) = &self.lw {
            writeln!(f, "Lw({carrier})={lw}")?;
        }

        if let Some((carrier, _, cn)) = &self.cn {
            writeln!(f, "Cn({carrier})={cn}")?;
        }

        Ok(())
    }
}

impl Difference {
    pub fn is_some(&self) -> bool {
        self.code.is_some()
            || self.code_j.is_some()
            || self.code_if.is_some()
            || self.phase.is_some()
            || self.phase_j.is_some()
            || self.phase_if.is_some()
            || self.cn.is_some()
            || self.lw.is_some()
    }

    pub fn with_code(mut self, value: (Carrier, f64)) -> Self {
        assert!(value.0.is_l1(), "internal error: this is not L1 frequency!");
        self.code = Some(value);
        self
    }

    pub fn with_code_j(mut self, value: (Carrier, f64)) -> Self {
        assert!(!value.0.is_l1(), "internal error: invalid frequency!");
        self.code_j = Some(value);
        self
    }

    pub fn with_code_if(mut self, value: (Carrier, f64)) -> Self {
        self.code_if = Some(value);
        self
    }

    pub fn with_phase(mut self, value: (Carrier, f64)) -> Self {
        assert!(value.0.is_l1(), "internal issue: this is not L1 frequency!");
        self.phase = Some(value);
        self
    }

    pub fn with_phase_j(mut self, value: (Carrier, f64)) -> Self {
        assert!(!value.0.is_l1(), "internal issue: this is L1 signal");
        self.phase_j = Some(value);
        self
    }

    pub fn with_phase_if(mut self, value: (Carrier, f64, f64)) -> Self {
        self.phase_if = Some(value);
        self
    }

    pub fn with_lw(mut self, value: (Carrier, f64, f64)) -> Self {
        self.lw = Some(value);
        self
    }

    pub fn with_cn(mut self, value: (Carrier, f64, f64)) -> Self {
        self.cn = Some(value);
        self
    }

    /// Obtain precise phase IF from purely differenced measurements.
    pub fn phase_if(&self, sv: SV) -> Option<f64> {
        let (c_1, phase_1) = self.phase?;
        let (c_j, phase_j) = self.phase_j?;

        let (lambda_1, lambda_j) = (
            SPEED_OF_LIGHT_M_S / c_1.frequency_hz(),
            SPEED_OF_LIGHT_M_S / c_j.frequency_hz(),
        );

        let (lambda_w, lambda_n) = (
            SPEED_OF_LIGHT_M_S / (c_1.frequency_hz() - c_j.frequency_hz()),
            SPEED_OF_LIGHT_M_S / (c_1.frequency_hz() + c_j.frequency_hz()),
        );

        let (_, _, phase_if) = self.phase_if?;

        let (_, _, phase_w) = self.lw?;
        let (_, _, code_n) = self.cn?;
        let nw = ((phase_w - code_n) / lambda_w).round();
        let n1 = ((phase_1 - phase_j - lambda_j * nw) / (lambda_1 - lambda_j)).round();
        let n2 = (n1 - nw).round();
        let nif = lambda_n * (n1 + lambda_w / lambda_j * nw);

        debug!("{sv} - Nw={nw}, N1={n1}, N2={n2}");
        debug!("{sv} - Nif={nif}");
        Some(phase_if - nif)
    }
}

/// [Differences] is a set to manage Single or Double [Difference]s.
#[derive(Default, Debug)]
pub(crate) struct Differences {
    /// [Difference] per [SV].
    pub inner: HashMap<SV, Difference>,
}

impl Differences {
    pub fn insert(&mut self, sv: SV, value: Difference) {
        self.inner.insert(sv, value);
    }

    /// Calculates the Double Difference (DD) from [Self] and rhs [Self].
    pub fn double_difference(&self, rhs: &Self) -> Self {
        let mut ret = Self::default();

        for (sat, difference) in self.inner.iter() {
            let mut diff = Difference::default();

            if let Some((_, rhs_diff)) = rhs.inner.iter().find(|(rhs_sat, _)| *rhs_sat == sat) {
                if let Some((c_1, p_1)) = difference.code {
                    if let Some((c_2, p_2)) = rhs_diff.code {
                        if c_1 == c_2 {
                            diff.code = Some((c_1, p_1 - p_2));
                        }
                    }
                }

                if let Some((c_1, p_1)) = difference.code_j {
                    if let Some((c_2, p_2)) = rhs_diff.code_j {
                        if c_1 == c_2 {
                            diff.code_j = Some((c_1, p_1 - p_2));
                        }
                    }
                }

                if let Some((c_1, p_1)) = difference.code_if {
                    if let Some((c_2, p_2)) = rhs_diff.code_if {
                        if c_1 == c_2 {
                            diff.code_if = Some((c_1, p_1 - p_2));
                        }
                    }
                }

                if let Some((c_1, p_1)) = difference.phase {
                    if let Some((c_2, p_2)) = rhs_diff.phase {
                        if c_1 == c_2 {
                            diff.phase = Some((c_1, p_1 - p_2));
                        }
                    }
                }

                if let Some((c_1, p_1)) = difference.phase_j {
                    if let Some((c_2, p_2)) = rhs_diff.phase_j {
                        if c_1 == c_2 {
                            diff.phase_j = Some((c_1, p_1 - p_2));
                        }
                    }
                }

                if let Some((c_1, lambda_1, p_1)) = difference.phase_if {
                    if let Some((c_2, _, p_2)) = rhs_diff.phase_if {
                        if c_1 == c_2 {
                            diff.phase_if = Some((c_1, lambda_1, p_1 - p_2));
                        }
                    }
                }

                if let Some((c_1, lambda_1, p_1)) = difference.lw {
                    if let Some((c_2, _, p_2)) = rhs_diff.lw {
                        if c_1 == c_2 {
                            diff.lw = Some((c_1, lambda_1, p_1 - p_2));
                        }
                    }
                }

                if let Some((c_1, lambda_1, p_1)) = difference.cn {
                    if let Some((c_2, _, p_2)) = rhs_diff.cn {
                        if c_1 == c_2 {
                            diff.cn = Some((c_1, lambda_1, p_1 - p_2));
                        }
                    }
                }
            }

            if diff.is_some() {
                ret.insert(*sat, diff);
            }
        }

        ret
    }

    /// Obtain [Difference] for this [SV] if it exists.
    pub fn difference(&self, sv: SV) -> Option<&Difference> {
        self.inner.get(&sv)
    }

    pub fn ndf(&self) -> usize {
        self.inner.len()
    }
}
