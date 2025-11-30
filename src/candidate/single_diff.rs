use log::debug;

use crate::{
    candidate::differences::Difference,
    // constants::SPEED_OF_LIGHT_M_S,
    prelude::Candidate,
};

impl Candidate {
    /// Runs the SD algorithm between [Self] and pivot [Self].
    pub(crate) fn single_difference(&self, pivot: &Self) -> Difference {
        let mut sd = Difference::default();

        if let Some((c_1, p_1)) = self.l1_pseudo_range() {
            if let Some((c_2, p_2)) = pivot.l1_pseudo_range() {
                if c_1 == c_2 {
                    sd = sd.with_code((c_1, p_1 - p_2));
                }
            }
        }

        if let Some((c_1, p_1)) = self.subsidary_pseudo_range() {
            if let Some((c_2, p_2)) = pivot.subsidary_pseudo_range() {
                if c_1 == c_2 {
                    sd = sd.with_code_j((c_1, p_1 - p_2));
                }
            }
        }

        if let Some(l_1) = self.code_if_combination() {
            if let Some(l_2) = pivot.code_if_combination() {
                if l_1.lhs == l_2.lhs && l_1.rhs == l_2.rhs {
                    sd = sd.with_code_if((l_1.rhs, l_1.value - l_2.value));
                }
            }
        }

        if let Some((c_1, l_1)) = self.l1_phase_range() {
            if let Some((c_2, l_2)) = pivot.l1_phase_range() {
                if c_1 == c_2 {
                    sd = sd.with_phase((c_1, l_1 - l_2));
                }
            }
        }

        if let Some((c_1, l_1)) = self.subsidary_phase_range() {
            if let Some((c_2, l_2)) = pivot.subsidary_phase_range() {
                if c_1 == c_2 {
                    sd = sd.with_phase_j((c_1, l_1 - l_2));
                }
            }
        }

        if let Some(l_1) = self.phase_if_combination() {
            if let Some(l_2) = pivot.phase_if_combination() {
                if l_1.lhs == l_2.lhs && l_1.rhs == l_2.rhs {
                    sd = sd.with_phase_if((l_1.rhs, l_1.lambda, l_1.value - l_2.value));
                }
            }
        }

        if let Some(c_1) = self.code_nl_combination() {
            if let Some(c_2) = pivot.code_nl_combination() {
                if c_1.lhs == c_2.lhs && c_1.rhs == c_2.rhs {
                    sd = sd.with_cn((c_1.rhs, c_1.lambda, c_1.value - c_2.value));
                }
            }
        }

        if let Some(l_1) = self.phase_wl_combination() {
            if let Some(l_2) = pivot.phase_wl_combination() {
                if l_1.lhs == l_2.lhs && l_1.rhs == l_2.rhs {
                    sd = sd.with_lw((l_1.rhs, l_1.lambda, l_1.value - l_2.value));
                }
            }
        }

        debug!("{}({}) - {}", self.epoch, self.sv, sd);
        sd
    }
}

#[cfg(test)]
mod test {

    use std::str::FromStr;

    use crate::{
        constants::SPEED_OF_LIGHT_M_S,
        prelude::{Carrier, Epoch},
        tests::{CandidatesBuilder, E01, E03, E05, init_logger},
    };

    #[test]
    fn single_difference_null() {
        let (f1, f2) = (Carrier::L1.frequency_hz(), Carrier::E5b.frequency_hz());
        let (f1pow, f2pow) = (f1.powi(2), f2.powi(2));
        let (lambda_1, lambda_j) = (SPEED_OF_LIGHT_M_S / f1, SPEED_OF_LIGHT_M_S / f2);

        let (lambda_n, lambda_w) = (
            SPEED_OF_LIGHT_M_S / (f1 + f2),
            SPEED_OF_LIGHT_M_S / (f1 - f2),
        );

        let f_if = f1 * f2 / (f1pow + f2pow).sqrt();
        let lambda_if = SPEED_OF_LIGHT_M_S / f_if;

        for t_str in [
            "2020-06-25T00:00:00 GPST",
            "2020-06-25T00:15:00 GPST",
            "2020-06-25T00:30:00 GPST",
            "2020-06-25T00:45:00 GPST",
            "2020-06-25T01:00:00 GPST",
        ] {
            let t = Epoch::from_str(t_str).unwrap();

            let pool = CandidatesBuilder::build_rover_at(t);

            assert!(
                !pool.is_empty(),
                "failed to propose any measurements at {t_str}"
            );

            for cd in pool.iter() {
                let single_diff = cd.single_difference(cd);

                assert_eq!(single_diff.code, Some((Carrier::L1, 0.0)));
                assert_eq!(single_diff.code_j, Some((Carrier::E5b, 0.0)));
                assert_eq!(single_diff.code_if, Some((Carrier::L1, 0.0)));

                assert_eq!(single_diff.phase, Some((Carrier::L1, 0.0)));
                assert_eq!(single_diff.phase_j, Some((Carrier::E5b, 0.0)));
                assert_eq!(single_diff.phase_if, Some((Carrier::L1, lambda_if, 0.0)));

                assert_eq!(single_diff.cn, Some((Carrier::L1, lambda_n, 0.0)));
                assert_eq!(single_diff.lw, Some((Carrier::L1, lambda_w, 0.0)));
            }
        }
    }

    #[test]
    fn rover_single_difference() {
        init_logger();

        let t_str = "2020-06-25T00:00:00 GPST";
        let t = Epoch::from_str(t_str).unwrap();

        let rover = CandidatesBuilder::build_rover_at(t);
        let e05 = CandidatesBuilder::build_rover_sv_at(E05, t);

        let (f1, fj) = (Carrier::L1.frequency_hz(), Carrier::E5b.frequency_hz());
        let (f1pow, fjpow) = (f1.powi(2), fj.powi(2));

        let (lambda_1, lambda_j) = (SPEED_OF_LIGHT_M_S / f1, SPEED_OF_LIGHT_M_S / fj);

        let (lambda_n, lambda_w) = (
            SPEED_OF_LIGHT_M_S / (f1 + fj),
            SPEED_OF_LIGHT_M_S / (f1 - fj),
        );

        let f_if = f1 * fj / (f1pow + fjpow).sqrt();
        let lambda_if = SPEED_OF_LIGHT_M_S / f_if;

        let mut e01_test_passed = false;
        let mut e03_test_passed = false;

        for rover in rover.iter() {
            let sd = rover.single_difference(&e05);

            if rover.sv == E01 {
                let p11 = 27616185.992;
                let p12 = 27616184.819;
                let p51 = 23730317.923;
                let p52 = 23730316.788;

                assert_eq!(
                    sd.code,
                    Some((Carrier::L1, p11 - p51)),
                    "{t_str}(E01/E05): invalid SD(P1-P1_ref)"
                );

                assert_eq!(
                    sd.code_j,
                    Some((Carrier::E5b, p12 - p52)),
                    "{t_str}(E01/E05): invalid SD(P2-P2_ref)"
                );

                let pc_1 = (f1pow * p11 - fjpow * p12) / (f1pow - fjpow);
                let pc_5 = (f1pow * p51 - fjpow * p52) / (f1pow - fjpow);

                assert_eq!(
                    sd.code_if,
                    Some((Carrier::L1, pc_1 - pc_5)),
                    "{t_str}(E01/E05): invalid SD(P1_if-P1_if_ref)"
                );

                let l11 = 145124050.106;
                let l12 = 108371872.760;
                let l51 = 124703702.220;
                let l52 = 93122915.921;

                assert_eq!(
                    sd.phase,
                    Some((Carrier::L1, l11 - l51)),
                    "{t_str}(E01/E05): invalid SD(L1-L1_ref)"
                );

                assert_eq!(
                    sd.phase_j,
                    Some((Carrier::E5b, l12 - l52)),
                    "{t_str}(E01/E05): invalid SD(L2-L2_ref)"
                );

                let lc_1 = (f1pow * l11 - fjpow * l12) / (f1pow - fjpow);
                let lc_5 = (f1pow * l51 - fjpow * l52) / (f1pow - fjpow);

                assert_eq!(
                    sd.phase_if,
                    Some((Carrier::L1, lambda_if, lc_1 - lc_5)),
                    "{t_str}(E01/E05): invalid SD(P1-P1_ref)"
                );

                let lw_1 = (f1 * l11 - fj * l12) / (f1 - fj);
                let lw_5 = (f1 * l51 - fj * l52) / (f1 - fj);

                assert_eq!(
                    sd.lw,
                    Some((Carrier::L1, lambda_w, lw_1 - lw_5)),
                    "{t_str}(E01/E05): invalid SD(Lw-Lw_ref)"
                );

                let cn_1 = (f1 * p11 + fj * p12) / (f1 + fj);
                let cn_5 = (f1 * p51 + fj * p52) / (f1 + fj);

                assert_eq!(
                    sd.cn,
                    Some((Carrier::L1, lambda_n, cn_1 - cn_5)),
                    "{t_str}(E01/E05): invalid SD(Cn-Cn_ref)"
                );

                e01_test_passed = true;
            } else if rover.sv == E03 {
                let p31 = 27055946.391;
                let p32 = 27055945.532;
                let p51 = 23730317.923;
                let p52 = 23730316.788;

                assert_eq!(
                    sd.code,
                    Some((Carrier::L1, p31 - p51)),
                    "{t_str}(E03/E05): invalid SD(L1-L1_ref)"
                );

                assert_eq!(
                    sd.code_j,
                    Some((Carrier::E5b, p32 - p52)),
                    "{t_str}(E03/E05): invalid SD(L2-L2_ref)"
                );

                let pc_3 = (f1pow * p31 - fjpow * p32) / (f1pow - fjpow);
                let pc_5 = (f1pow * p51 - fjpow * p52) / (f1pow - fjpow);

                assert_eq!(
                    sd.code_if,
                    Some((Carrier::L1, pc_3 - pc_5)),
                    "{t_str}(E03/E05): invalid SD(P1-P1_ref)"
                );

                let l31 = 142179967.778;
                let l32 = 106173364.686;
                let l51 = 124703702.220;
                let l52 = 93122915.921;

                assert_eq!(
                    sd.phase,
                    Some((Carrier::L1, l31 - l51)),
                    "{t_str}(E03/E05): invalid SD(L1-L1_ref)"
                );

                assert_eq!(
                    sd.phase_j,
                    Some((Carrier::E5b, l32 - l52)),
                    "{t_str}(E03/E05): invalid SD(L2-L2_ref)"
                );

                let lc_3 = (f1pow * l31 - fjpow * l32) / (f1pow - fjpow);
                let lc_5 = (f1pow * l51 - fjpow * l52) / (f1pow - fjpow);

                assert_eq!(
                    sd.phase_if,
                    Some((Carrier::L1, lambda_if, lc_3 - lc_5)),
                    "{t_str}(E03/E05): invalid SD(Lif-L2if_ref)"
                );

                let lw_3 = (f1 * l31 - fj * l32) / (f1 - fj);
                let lw_5 = (f1 * l51 - fj * l52) / (f1 - fj);

                assert_eq!(
                    sd.lw,
                    Some((Carrier::L1, lambda_w, lw_3 - lw_5)),
                    "{t_str}(E03/E05): invalid SD(Lw-Lw_ref)"
                );

                let cn_3 = (f1 * p31 + fj * p32) / (f1 + fj);
                let cn_5 = (f1 * p51 + fj * p52) / (f1 + fj);

                assert_eq!(
                    sd.cn,
                    Some((Carrier::L1, lambda_n, cn_3 - cn_5)),
                    "{t_str}(E03/E05): invalid SD(Cn-Cn_ref)"
                );

                e03_test_passed = true;
            }
        }

        assert!(e01_test_passed);
        assert!(e03_test_passed);
    }
}
