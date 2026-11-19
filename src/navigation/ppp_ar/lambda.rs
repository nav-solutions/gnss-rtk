//! LAMBDA / MLAMBDA integer least-squares ambiguity resolution.
//!
//! This is a faithful port of the reference implementation found in RTKLIB
//! (`lambda.c`, T. Takasu), itself based on:
//!
//! - [1] P.J.G. Teunissen, "The least-squares ambiguity decorrelation adjustment:
//!   a method for fast GPS ambiguity estimation", J. Geodesy, 1995.
//! - [2] X.-W. Chang, X. Yang, T. Zhou, "MLAMBDA: a modified LAMBDA method for
//!   integer least-squares estimation", J. Geodesy, 2005.
//!
//! Conventions (identical to RTKLIB):
//! - `Q = Lᵀ · diag(D) · L`, `L` unit lower triangular.
//! - Decorrelation `Z` (unimodular) so that `z = Zᵀ a`, `Qz = Zᵀ Q Z`.
//! - The `m` best integer candidates are returned as columns of `F`, sorted
//!   by increasing squared residual norm `s` (`s[0]` is the best).
use nalgebra::{DMatrix, DVector};

use crate::prelude::Error;

use log::debug;

pub struct LambdaAR {}

impl LambdaAR {
    /// Maximal number of search loops
    pub const MAX_SEARCH: usize = 10_000;

    fn signum(value: f64) -> f64 {
        if value <= 0.0 { -1.0 } else { 1.0 }
    }

    fn round(value: f64) -> f64 {
        (value + 0.5).floor()
    }

    /// LD factorization: `Q = Lᵀ · diag(D) · L`.
    pub(crate) fn ld_factorization(
        q_mat: &DMatrix<f64>,
    ) -> Result<(DMatrix<f64>, DVector<f64>), Error> {
        let n = q_mat.nrows();

        let mut a = q_mat.clone();
        let mut l_mat = DMatrix::<f64>::zeros(n, n);
        let mut d_vec = DVector::<f64>::zeros(n);

        for i in (0..n).rev() {
            d_vec[i] = a[(i, i)];

            if d_vec[i] <= 0.0 {
                return Err(Error::AmbiguityFactorization);
            }

            let sq = d_vec[i].sqrt();

            for j in 0..=i {
                l_mat[(i, j)] = a[(i, j)] / sq;
            }

            for j in 0..i {
                for k in 0..=j {
                    a[(j, k)] -= l_mat[(i, k)] * l_mat[(i, j)];
                }
            }

            let lii = l_mat[(i, i)];

            for j in 0..=i {
                l_mat[(i, j)] /= lii;
            }
        }

        Ok((l_mat, d_vec))
    }

    /// Integer Gauss transformation
    pub(crate) fn gauss_transform(
        i: usize,
        j: usize,
        ndf: usize,
        l_mat: &mut DMatrix<f64>,
        z_mat: &mut DMatrix<f64>,
    ) {
        let mu = Self::round(l_mat[(i, j)]);

        if mu != 0.0 {
            for k in i..ndf {
                l_mat[(k, j)] -= mu * l_mat[(k, i)];
            }

            for k in 0..ndf {
                z_mat[(k, j)] -= mu * z_mat[(k, i)];
            }
        }
    }

    /// Permutation of (j, j+1)
    fn permutations(
        ndf: usize,
        l_mat: &mut DMatrix<f64>,
        d_vec: &mut DVector<f64>,
        j: usize,
        delta: f64,
        z_mat: &mut DMatrix<f64>,
    ) {
        let eta = d_vec[j] / delta;
        let lambda = d_vec[j + 1] * l_mat[(j + 1, j)] / delta;

        d_vec[j] = eta * d_vec[j + 1];
        d_vec[j + 1] = delta;

        for k in 0..j {
            let a0 = l_mat[(j, k)];
            let a1 = l_mat[(j + 1, k)];

            l_mat[(j, k)] = -l_mat[(j + 1, j)] * a0 + a1;
            l_mat[(j + 1, k)] = eta * a0 + lambda * a1;
        }

        l_mat[(j + 1, j)] = lambda;

        for k in j + 2..ndf {
            l_mat.swap((k, j), (k, j + 1));
        }

        for k in 0..ndf {
            z_mat.swap((k, j), (k, j + 1));
        }
    }

    /// LAMBDA reduction (decorrelation + sorting).
    fn reduction(
        ndf: usize,
        l_mat: &mut DMatrix<f64>,
        d_vec: &mut DVector<f64>,
        z_mat: &mut DMatrix<f64>,
    ) {
        if ndf < 2 {
            return;
        }

        // signed indexes: loop runs down to j = 0 included
        let mut j = ndf as isize - 2;
        let mut k = ndf as isize - 2;

        while j >= 0 {
            let ju = j as usize;

            if j <= k {
                for i in ju + 1..ndf {
                    Self::gauss_transform(i, ju, ndf, l_mat, z_mat);
                }
            }

            let delta = d_vec[ju] + l_mat[(ju + 1, ju)] * l_mat[(ju + 1, ju)] * d_vec[ju + 1];

            if delta + 1E-6 < d_vec[ju + 1] {
                // compared considering numerical error
                Self::permutations(ndf, l_mat, d_vec, ju, delta, z_mat);
                k = j;
                j = ndf as isize - 2;
            } else {
                j -= 1;
            }
        }
    }

    /// MLAMBDA search: `m` best candidates, returned in `zn` (n x m), and their
    /// squared residual norms in `s` (sorted, best first).
    fn search(
        ndf: usize,
        nfixed: usize,
        l_mat: &DMatrix<f64>,
        d_vec: &DVector<f64>,
        zs_vec: &DVector<f64>,
        zn_mat: &mut DMatrix<f64>,
        s_vec: &mut DVector<f64>,
    ) -> Result<(), Error> {
        let n = ndf;
        let m = nfixed;

        let mut nn = 0usize;
        let mut imax = 0usize;
        let mut maxdist = 1E99_f64;

        let mut s_mat = DMatrix::<f64>::zeros(n, n);
        let mut dist = DVector::<f64>::zeros(n);
        let mut zb = DVector::<f64>::zeros(n);
        let mut z = DVector::<f64>::zeros(n);
        let mut step = DVector::<f64>::zeros(n);

        let mut k = n - 1;

        dist[k] = 0.0;
        zb[k] = zs_vec[k];
        z[k] = Self::round(zb[k]);

        let mut y = zb[k] - z[k];
        step[k] = Self::signum(y);

        let mut c = 0;

        while c < Self::MAX_SEARCH {
            let newdist = dist[k] + y * y / d_vec[k];

            if newdist < maxdist {
                if k != 0 {
                    // Case 1: move down
                    k -= 1;
                    dist[k] = newdist;

                    for i in 0..=k {
                        s_mat[(k, i)] =
                            s_mat[(k + 1, i)] + (z[k + 1] - zb[k + 1]) * l_mat[(k + 1, i)];
                    }

                    zb[k] = zs_vec[k] + s_mat[(k, k)];
                    z[k] = Self::round(zb[k]);
                    y = zb[k] - z[k];
                    step[k] = Self::signum(y);
                } else {
                    // Case 2: store the candidate and try next valid integer
                    if nn < m {
                        if nn == 0 || newdist > s_vec[imax] {
                            imax = nn;
                        }

                        for i in 0..n {
                            zn_mat[(i, nn)] = z[i];
                        }

                        s_vec[nn] = newdist;
                        nn += 1;
                    } else {
                        if newdist < s_vec[imax] {
                            for i in 0..n {
                                zn_mat[(i, imax)] = z[i];
                            }

                            s_vec[imax] = newdist;

                            imax = 0;
                            for i in 0..m {
                                if s_vec[imax] < s_vec[i] {
                                    imax = i;
                                }
                            }
                        }

                        maxdist = s_vec[imax];
                    }

                    z[0] += step[0]; // next valid integer
                    y = zb[0] - z[0];
                    step[0] = -step[0] - Self::signum(step[0]);
                }
            } else {
                // Case 3: exit or move up
                if k == n - 1 {
                    break;
                }

                k += 1;
                z[k] += step[k];
                y = zb[k] - z[k];
                step[k] = -step[k] - Self::signum(step[k]);
            }

            c += 1;
        }

        // sort by s
        for i in 0..m.saturating_sub(1) {
            for j in i + 1..m {
                if s_vec[i] < s_vec[j] {
                    continue;
                }

                s_vec.swap_rows(i, j);

                for k in 0..n {
                    zn_mat.swap((k, i), (k, j));
                }
            }
        }

        if c >= Self::MAX_SEARCH {
            return Err(Error::UnresolvedAmbiguity);
        }

        if nn < m {
            // less candidates than requested were found
            return Err(Error::UnresolvedAmbiguity);
        }

        Ok(())
    }

    /// Runs the LAMBDA / MLAMBDA integer least-squares estimation.
    ///
    /// ## Input
    /// - ndf: number of float ambiguities (n)
    /// - nfixed: number of integer candidates to return (m)
    /// - x_vec: float ambiguities (n x 1)
    /// - q_mat: float ambiguities covariance (n x n), symmetric positive definite
    ///
    /// ## Output
    /// - F (n x m): integer candidates, as columns, best first.
    /// - s (m): squared residual norms `(x - F_i)ᵀ Q⁻¹ (x - F_i)`, sorted (best first).
    pub fn run(
        ndf: usize,
        nfixed: usize,
        x_vec: &DMatrix<f64>,
        q_mat: &DMatrix<f64>,
    ) -> Result<(DMatrix<f64>, DVector<f64>), Error> {
        let (x_rows, x_cols) = (x_vec.nrows(), x_vec.ncols());
        let (q_rows, q_cols) = (q_mat.nrows(), q_mat.ncols());

        if ndf == 0 || nfixed == 0 {
            return Err(Error::MatrixMinimalDimension);
        }

        if x_cols != 1 || x_rows != ndf || q_rows != ndf || q_cols != ndf {
            return Err(Error::MatrixDimension);
        }

        let a = DVector::<f64>::from_column_slice(x_vec.as_slice());

        debug!("(ppp) lambda - ndf={ndf} - X={x_vec} Q={q_mat}");

        let (mut l_mat, mut d_vec) = Self::ld_factorization(q_mat)?;
        let mut z_mat = DMatrix::<f64>::identity(ndf, ndf);

        Self::reduction(ndf, &mut l_mat, &mut d_vec, &mut z_mat);

        // z = Zᵀ a
        let zs_vec = z_mat.transpose() * &a;

        let mut e_mat = DMatrix::<f64>::zeros(ndf, nfixed);
        let mut s_vec = DVector::<f64>::zeros(nfixed);

        Self::search(ndf, nfixed, &l_mat, &d_vec, &zs_vec, &mut e_mat, &mut s_vec)?;

        // F = Zᵀ \ E
        let zt_inv = z_mat
            .transpose()
            .try_inverse()
            .ok_or(Error::AmbiguityInverse)?;

        let mut f_mat = zt_inv * e_mat;

        // Z is unimodular: the result is integer, remove numerical noise.
        f_mat.iter_mut().for_each(|v| *v = v.round());

        debug!("(ppp) lambda - F={f_mat} S={s_vec}");

        Ok((f_mat, s_vec))
    }

    /// Ratio test (second best / best squared norms).
    /// Returns None when less than two candidates are available or best norm is null.
    pub fn ratio(s_vec: &DVector<f64>) -> Option<f64> {
        if s_vec.nrows() < 2 {
            return None;
        }

        if s_vec[0] <= 0.0 {
            // perfect match: infinitely confident
            return Some(f64::INFINITY);
        }

        Some(s_vec[1] / s_vec[0])
    }
}

#[cfg(test)]
mod test {

    use super::LambdaAR;
    use crate::tests::init_logger;

    use nalgebra::{DMatrix, DVector, DimName, U1, U6, U10};

    /// Weighted squared norm (a - z)ᵀ Q⁻¹ (a - z)
    fn sq_norm(a: &DVector<f64>, z: &DVector<f64>, q_inv: &DMatrix<f64>) -> f64 {
        let d = a - z;
        (d.transpose() * q_inv * &d)[(0, 0)]
    }

    /// Brute force integer least squares, used as an independent oracle.
    /// Enumerates all integer vectors within +/- `radius` of the rounded float solution
    /// and returns the 2 best candidates and their norms.
    fn brute_force_ils(
        a: &DVector<f64>,
        q: &DMatrix<f64>,
        radius: i64,
    ) -> Vec<(DVector<f64>, f64)> {
        let n = a.nrows();
        let q_inv = q.clone().try_inverse().unwrap();
        let q_inv: Vec<f64> = (0..n * n).map(|k| q_inv[(k / n, k % n)]).collect();
        let center: Vec<i64> = a.iter().map(|v| v.round() as i64).collect();

        let width = (2 * radius + 1) as usize;
        let total = width.pow(n as u32);

        let mut best: Vec<(Vec<f64>, f64)> = Vec::with_capacity(3);
        let mut z = vec![0.0; n];
        let mut d = vec![0.0; n];

        for idx in 0..total {
            let mut rem = idx;

            for i in 0..n {
                let off = (rem % width) as i64 - radius;
                rem /= width;
                z[i] = (center[i] + off) as f64;
                d[i] = a[i] - z[i];
            }

            let mut norm = 0.0;
            for i in 0..n {
                let mut acc = 0.0;
                for j in 0..n {
                    acc += q_inv[i * n + j] * d[j];
                }
                norm += d[i] * acc;
            }

            if best.len() < 2 || norm < best[1].1 {
                best.push((z.clone(), norm));
                best.sort_by(|x, y| x.1.partial_cmp(&y.1).unwrap());
                best.truncate(2);
            }
        }

        best.into_iter()
            .map(|(z, s)| (DVector::from_vec(z), s))
            .collect()
    }

    /// Deterministic pseudo random generator (LCG), no external dependency.
    struct Lcg(u64);

    impl Lcg {
        fn next_f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((self.0 >> 11) as f64) / ((1u64 << 53) as f64)
        }

        fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
            lo + (hi - lo) * self.next_f64()
        }
    }

    /// Builds a random, highly correlated SPD matrix, typical of GNSS float ambiguities:
    /// Q = s² (A Aᵀ + eps I) with A having a dominant common component.
    fn random_gnss_like_q(rng: &mut Lcg, n: usize, scale: f64) -> DMatrix<f64> {
        let mut a = DMatrix::<f64>::zeros(n, n + 3);

        for i in 0..n {
            // common "geometry" components → strong correlations
            for j in 0..3 {
                a[(i, j)] = rng.uniform(0.5, 1.5);
            }
            for j in 3..n + 3 {
                a[(i, j)] = rng.uniform(-0.1, 0.1);
            }
        }

        let mut q = &a * a.transpose();

        for i in 0..n {
            q[(i, i)] += 1.0E-3;
        }

        q * scale
    }

    #[test]
    fn gauss_transform() {
        let mut l_mat = DMatrix::<f64>::identity(U6::DIM, U6::DIM);
        let mut z_mat = l_mat.clone();

        for i in 0..U6::DIM {
            for j in 0..U6::DIM {
                LambdaAR::gauss_transform(i, j, U6::DIM, &mut l_mat, &mut z_mat);
            }
        }
    }

    #[test]
    fn ld_factorization_reconstructs_q() {
        let mut rng = Lcg(42);

        for n in [1, 2, 3, 5, 8, 12] {
            let q = random_gnss_like_q(&mut rng, n, 0.3);
            let (l, d) = LambdaAR::ld_factorization(&q).unwrap();

            // L unit lower triangular
            for i in 0..n {
                assert!((l[(i, i)] - 1.0).abs() < 1E-12);
                for j in i + 1..n {
                    assert_eq!(l[(i, j)], 0.0);
                }
            }

            let rebuilt = l.transpose() * DMatrix::from_diagonal(&d) * &l;

            for i in 0..n {
                for j in 0..n {
                    assert!(
                        (rebuilt[(i, j)] - q[(i, j)]).abs() < 1E-9,
                        "n={n} LtDL mismatch at ({i},{j})"
                    );
                }
            }
        }
    }

    #[test]
    fn ld_factorization_rejects_non_positive_definite() {
        let q = DMatrix::<f64>::from_row_slice(2, 2, &[1.0, 2.0, 2.0, 1.0]);
        assert!(LambdaAR::ld_factorization(&q).is_err());
    }

    #[test]
    fn invalid_dimensions() {
        let x = DMatrix::<f64>::zeros(3, 1);
        let q = DMatrix::<f64>::identity(2, 2);
        assert!(LambdaAR::run(3, 2, &x, &q).is_err());
        assert!(LambdaAR::run(0, 2, &x, &q).is_err());
        assert!(LambdaAR::run(3, 0, &x, &DMatrix::identity(3, 3)).is_err());
    }

    #[test]
    fn single_ambiguity_is_rounding() {
        for (a, expected_1, expected_2) in [
            (1.2, 1.0, 2.0),
            (-3.7, -4.0, -3.0),
            (10.49, 10.0, 11.0),
            (-0.51, -1.0, 0.0),
        ] {
            let x = DMatrix::<f64>::from_row_slice(1, 1, &[a]);
            let q = DMatrix::<f64>::from_row_slice(1, 1, &[0.04]);
            let (f, s) = LambdaAR::run(1, 2, &x, &q).unwrap();
            assert_eq!(f[(0, 0)], expected_1, "a={a}");
            assert_eq!(f[(0, 1)], expected_2, "a={a}");
            assert!((s[0] - (a - expected_1).powi(2) / 0.04).abs() < 1E-9);
            assert!(s[0] <= s[1]);
        }
    }

    #[test]
    fn integer_input_is_fixed_exactly() {
        let mut rng = Lcg(7);

        for n in [2, 4, 7, 10] {
            let q = random_gnss_like_q(&mut rng, n, 0.05);
            let truth: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0E6, 1.0E6).round()).collect();
            let x = DMatrix::<f64>::from_column_slice(n, 1, &truth);

            let (f, s) = LambdaAR::run(n, 2, &x, &q).unwrap();

            for i in 0..n {
                assert_eq!(f[(i, 0)], truth[i], "n={n} i={i}");
            }

            assert!(s[0].abs() < 1E-9, "best norm should be null");
            assert!(s[1] > 0.0);
            assert_eq!(LambdaAR::ratio(&s), Some(f64::INFINITY));
        }
    }

    /// RTKLIB unit test #1 (t_lambda.c, utest1)
    #[test]
    fn mlambda_rtklib_utest1() {
        init_logger();

        let x = DMatrix::<f64>::from_row_slice(
            U6::DIM,
            U1::DIM,
            &[
                1585184.171,
                -6716599.430,
                3915742.905,
                7627233.455,
                9565990.879,
                989457273.200,
            ],
        );

        let q = DMatrix::<f64>::from_row_slice(
            U6::DIM,
            U6::DIM,
            &[
                0.227134, 0.112202, 0.112202, 0.112202, 0.112202, 0.103473, 0.112202, 0.227134,
                0.112202, 0.112202, 0.112202, 0.103473, 0.112202, 0.112202, 0.227134, 0.112202,
                0.112202, 0.103473, 0.112202, 0.112202, 0.112202, 0.227134, 0.112202, 0.103473,
                0.112202, 0.112202, 0.112202, 0.112202, 0.227134, 0.103473, 0.103473, 0.103473,
                0.103473, 0.103473, 0.103473, 0.434339,
            ],
        );

        // expected candidates, as columns (best first)
        let f_expected = DMatrix::<f64>::from_row_slice(
            U6::DIM,
            2,
            &[
                1585184.0,
                1585184.0,
                -6716599.0,
                -6716600.0,
                3915743.0,
                3915743.0,
                7627234.0,
                7627233.0,
                9565991.0,
                9565991.0,
                989457273.0,
                989457273.0,
            ],
        );

        let s_expected = DVector::<f64>::from_row_slice(&[3.507984, 3.708456]);

        let (f, s) = LambdaAR::run(U6::DIM, 2, &x, &q).unwrap_or_else(|e| {
            panic!("mlambda search failed with {e}");
        });

        assert_eq!(f, f_expected, "invalid candidates");

        for i in 0..2 {
            assert!(
                (s[i] - s_expected[i]).abs() < 1E-4,
                "s[{i}]={} expected {}",
                s[i],
                s_expected[i]
            );
        }

        let ratio = LambdaAR::ratio(&s).unwrap();
        assert!((ratio - 3.708456 / 3.507984).abs() < 1E-4);
    }

    /// RTKLIB unit test #2 (t_lambda.c, utest2)
    #[test]
    fn mlambda_rtklib_utest2() {
        init_logger();

        let a = DMatrix::<f64>::from_row_slice(
            U10::DIM,
            U1::DIM,
            &[
                -13324172.755747,
                -10668894.713608,
                -7157225.010770,
                -6149367.974367,
                -7454133.571066,
                -5969200.494550,
                8336734.058423,
                6186974.084502,
                -17549093.883655,
                -13970158.922370,
            ],
        );

        let q = DMatrix::<f64>::from_row_slice(
            U10::DIM,
            U10::DIM,
            &[
                0.446320, 0.223160, 0.223160, 0.223160, 0.223160, 0.572775, 0.286388, 0.286388,
                0.286388, 0.286388, 0.223160, 0.446320, 0.223160, 0.223160, 0.223160, 0.286388,
                0.572775, 0.286388, 0.286388, 0.286388, 0.223160, 0.223160, 0.446320, 0.223160,
                0.223160, 0.286388, 0.286388, 0.572775, 0.286388, 0.286388, 0.223160, 0.223160,
                0.223160, 0.446320, 0.223160, 0.286388, 0.286388, 0.286388, 0.572775, 0.286388,
                0.223160, 0.223160, 0.223160, 0.223160, 0.446320, 0.286388, 0.286388, 0.286388,
                0.286388, 0.572775, 0.572775, 0.286388, 0.286388, 0.286388, 0.286388, 0.735063,
                0.367531, 0.367531, 0.367531, 0.367531, 0.286388, 0.572775, 0.286388, 0.286388,
                0.286388, 0.367531, 0.735063, 0.367531, 0.367531, 0.367531, 0.286388, 0.286388,
                0.572775, 0.286388, 0.286388, 0.367531, 0.367531, 0.735063, 0.367531, 0.367531,
                0.286388, 0.286388, 0.286388, 0.572775, 0.286388, 0.367531, 0.367531, 0.367531,
                0.735063, 0.367531, 0.286388, 0.286388, 0.286388, 0.286388, 0.572775, 0.367531,
                0.367531, 0.367531, 0.367531, 0.735063,
            ],
        );

        let f_expected = DMatrix::<f64>::from_row_slice(
            U10::DIM,
            2,
            &[
                -13324188.0,
                -13324188.0,
                -10668901.0,
                -10668908.0,
                -7157236.0,
                -7157236.0,
                -6149379.0,
                -6149379.0,
                -7454143.0,
                -7454143.0,
                -5969220.0,
                -5969220.0,
                8336726.0,
                8336717.0,
                6186960.0,
                6186960.0,
                -17549108.0,
                -17549108.0,
                -13970171.0,
                -13970171.0,
            ],
        );

        let s_expected = DVector::<f64>::from_row_slice(&[1506.435789, 1612.811795]);

        let (f, s) = LambdaAR::run(U10::DIM, 2, &a, &q).unwrap_or_else(|e| {
            panic!("mlambda search failed with {e}");
        });

        assert_eq!(f, f_expected, "invalid candidates");

        for i in 0..2 {
            assert!(
                (s[i] - s_expected[i]).abs() < 1E-3,
                "s[{i}]={} expected {}",
                s[i],
                s_expected[i]
            );
        }
    }

    /// Classical 3D example (Teunissen / de Jonge & Tiberius LAMBDA tutorial),
    /// checked against brute force.
    #[test]
    fn lambda_textbook_3d() {
        let a = DMatrix::<f64>::from_row_slice(3, 1, &[5.45, 3.10, 2.97]);

        let q = DMatrix::<f64>::from_row_slice(
            3,
            3,
            &[
                6.290, 5.978, 0.544, 5.978, 6.292, 2.340, 0.544, 2.340, 6.288,
            ],
        );

        let (f, s) = LambdaAR::run(3, 2, &a, &q).unwrap();

        let oracle = brute_force_ils(&DVector::from_column_slice(a.as_slice()), &q, 6);

        for k in 0..2 {
            for i in 0..3 {
                assert_eq!(f[(i, k)], oracle[k].0[i], "candidate #{k}");
            }
            assert!((s[k] - oracle[k].1).abs() < 1E-9);
        }

        // decorrelation must beat plain rounding here
        let rounded = DVector::<f64>::from_row_slice(&[5.0, 3.0, 3.0]);
        let q_inv = q.clone().try_inverse().unwrap();
        assert!(s[0] <= sq_norm(&DVector::from_column_slice(a.as_slice()), &rounded, &q_inv));
    }

    /// Randomized comparison against brute force integer least squares.
    #[test]
    fn mlambda_vs_brute_force() {
        let mut rng = Lcg(2026);

        for trial in 0..60 {
            let n = 2 + (trial % 3); // 2..=4
            let scale = rng.uniform(0.01, 0.5);
            let q = random_gnss_like_q(&mut rng, n, scale);

            let a: Vec<f64> = (0..n).map(|_| rng.uniform(-50.0, 50.0)).collect();
            let a_vec = DVector::from_column_slice(&a);
            let x = DMatrix::<f64>::from_column_slice(n, 1, &a);

            let (f, s) = LambdaAR::run(n, 2, &x, &q).unwrap_or_else(|e| {
                panic!("trial {trial}: lambda failed: {e}");
            });

            // search radius large enough with respect to the float uncertainty
            let radius =
                (3.0 * (0..n).map(|i| q[(i, i)].sqrt()).fold(0.0, f64::max)).ceil() as i64 + 2;
            let radius = radius.min(5);

            let oracle = brute_force_ils(&a_vec, &q, radius);

            assert!(
                (s[0] - oracle[0].1).abs() < 1E-6 * (1.0 + oracle[0].1),
                "trial {trial} (n={n}): best norm {} vs oracle {}",
                s[0],
                oracle[0].1
            );

            assert!(
                (s[1] - oracle[1].1).abs() < 1E-6 * (1.0 + oracle[1].1),
                "trial {trial} (n={n}): 2nd norm {} vs oracle {}",
                s[1],
                oracle[1].1
            );

            for i in 0..n {
                assert_eq!(
                    f[(i, 0)],
                    oracle[0].0[i],
                    "trial {trial} (n={n}) best candidate"
                );
            }
        }
    }

    /// Candidates are sorted and norms are consistent with the returned integers.
    #[test]
    fn candidates_sorted_and_consistent() {
        let mut rng = Lcg(99);

        for n in [3, 6, 9, 14] {
            let q = random_gnss_like_q(&mut rng, n, 0.2);
            let q_inv = q.clone().try_inverse().unwrap();
            let a: Vec<f64> = (0..n).map(|_| rng.uniform(-1.0E4, 1.0E4)).collect();
            let a_vec = DVector::from_column_slice(&a);
            let x = DMatrix::<f64>::from_column_slice(n, 1, &a);

            let m = 5;
            let (f, s) = LambdaAR::run(n, m, &x, &q).unwrap();

            for k in 0..m {
                let z = f.column(k).into_owned();
                let norm = sq_norm(&a_vec, &z, &q_inv);
                assert!((norm - s[k]).abs() < 1E-6 * (1.0 + norm), "n={n} k={k}");
                if k > 0 {
                    assert!(s[k - 1] <= s[k], "n={n}: not sorted");
                    assert_ne!(f.column(k - 1), f.column(k), "n={n}: duplicate candidates");
                }
            }
        }
    }

    /// Ambiguities known to be integers + small noise: LAMBDA must recover the truth
    /// with a large ratio, even with strong correlations where rounding fails.
    #[test]
    fn mlambda_recovers_truth_on_correlated_noise() {
        let mut rng = Lcg(31415);
        let mut success = 0;
        let trials = 50;

        for _ in 0..trials {
            let n = 8;
            let q = random_gnss_like_q(&mut rng, n, 0.002);

            // correlated noise sample: x = truth + Lc · w, with Q = Lc Lcᵀ
            let chol = q.clone().cholesky().unwrap();
            let w = DVector::<f64>::from_iterator(
                n,
                (0..n).map(|_| {
                    // Box-Muller
                    let u1 = rng.next_f64().max(1E-12);
                    let u2 = rng.next_f64();
                    (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
                }),
            );

            let truth =
                DVector::<f64>::from_iterator(n, (0..n).map(|_| rng.uniform(-1E5, 1E5).round()));
            let noisy = &truth + chol.l() * w;

            let x = DMatrix::<f64>::from_column_slice(n, 1, noisy.as_slice());
            let (f, s) = LambdaAR::run(n, 2, &x, &q).unwrap();

            if f.column(0) == truth.column(0) {
                success += 1;
                assert!(LambdaAR::ratio(&s).unwrap() >= 1.0);
            }
        }

        assert!(
            success >= trials * 9 / 10,
            "success rate too low: {success}/{trials}"
        );
    }
}
