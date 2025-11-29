//! Brancroft solver
use crate::{
    constants::{EARTH_EQUATORIAL_RADIUS_KM, SPEED_OF_LIGHT_M_S},
    error::Error,
    prelude::Candidate,
};
use log::error;

use nalgebra::{Matrix4, Vector4};

pub struct Bancroft {
    a: Vector4<f64>,
    b: Matrix4<f64>,
    m: Matrix4<f64>,
    ones: Vector4<f64>,
}

fn lorentz_4_4(a: Vector4<f64>, b: Vector4<f64>, m: &Matrix4<f64>) -> f64 {
    let scalar = a.transpose() * m * b;
    scalar[(0, 0)]
}

impl Bancroft {
    /// Generates the m [Matrix4]
    fn m_matrix() -> Matrix4<f64> {
        let mut m = Matrix4::<f64>::identity();
        m[(3, 3)] = -1.0;
        m
    }

    /// Generates a [Vector4] with 1's
    fn one_vector() -> Vector4<f64> {
        Vector4::<f64>::new(1.0_f64, 1.0_f64, 1.0_f64, 1.0_f64)
    }

    /// Builds new Bancroft solver
    pub fn new(cd: &[Candidate]) -> Result<Self, Error> {
        let m = Self::m_matrix();
        let mut a = Vector4::<f64>::default();
        let mut b = Matrix4::<f64>::default();

        if cd.len() < 4 {
            return Err(Error::NotEnoughInitializationCandidates);
        }

        let mut j = 0;

        for i in 0..cd.len() {
            if let Some(orbit) = cd[i].orbit {
                let state = orbit.to_cartesian_pos_vel() * 1.0E3;
                let (x_i, y_i, z_i) = (state[0], state[1], state[2]);

                if let Some((_, r_i)) = cd[i].best_snr_range_m() {
                    let dt_i = cd[i].clock_corr.duration.to_seconds();
                    let tgd_i = cd[i].tgd.to_seconds();
                    let pr_i = r_i + dt_i * SPEED_OF_LIGHT_M_S - tgd_i;

                    b[(j, 0)] = x_i;
                    b[(j, 1)] = y_i;
                    b[(j, 2)] = z_i;
                    b[(j, 3)] = pr_i;
                    a[j] = 0.5 * (x_i.powi(2) + y_i.powi(2) + z_i.powi(2) - pr_i.powi(2));

                    j += 1;

                    if j == 4 {
                        break;
                    }
                }
            } else {
                error!(
                    "{}({}) bancroft unresolved orbital state",
                    cd[i].epoch, cd[i].sv
                );
            }
        }
        if j != 4 {
            Err(Error::BancroftError)
        } else {
            Ok(Self {
                a,
                b,
                m,
                ones: Self::one_vector(),
            })
        }
    }

    /// [Bancroft] resolution
    pub fn resolve(&self) -> Result<Vector4<f64>, Error> {
        let r_e = EARTH_EQUATORIAL_RADIUS_KM * 1.0E3;

        let b_inv = self.b.try_inverse().ok_or(Error::MatrixInversion)?;

        let b_1 = b_inv * self.ones;
        let b_a = b_inv * self.a;

        let a = lorentz_4_4(b_1, b_1, &self.m);
        let b = 2.0 * (lorentz_4_4(b_1, b_a, &self.m) - 1.0);

        let c = lorentz_4_4(b_a, b_a, &self.m);

        let delta = b.powi(2) - 4.0 * a * c;

        if delta > 0.0 {
            let delta_sqrt = delta.sqrt();
            let x = ((-b + delta_sqrt) / 2.0 / a, (-b - delta_sqrt) / 2.0 / a);

            let solutions = (
                self.m * b_inv * (x.0 * self.ones + self.a),
                self.m * b_inv * (x.1 * self.ones + self.a),
            );

            let rho = (
                (solutions.0[0].powi(2) + solutions.0[1].powi(2) + solutions.0[2].powi(2)).sqrt(),
                (solutions.1[0].powi(2) + solutions.1[1].powi(2) + solutions.1[2].powi(2)).sqrt(),
            );

            let err = ((rho.0 - r_e).abs(), (rho.1 - r_e).abs());

            if err.0 < err.1 {
                Ok(solutions.0)
            } else {
                Ok(solutions.1)
            }
        } else if delta < 0.0 {
            Err(Error::BancroftImaginarySolution)
        } else {
            let x = -b / a / 2.0;
            Ok(self.m * b_inv * (x * self.ones + self.a))
        }
    }
}

#[cfg(test)]
mod test {
    use super::{Bancroft, lorentz_4_4};
    use nalgebra::Vector4;

    #[test]
    fn lorentz_4_4_product() {
        let a = Vector4::<f64>::new(1.0, 2.0, 3.0, 4.0);
        let b = Vector4::<f64>::new(5.0, 6.0, 7.0, 8.0);
        let m = Bancroft::m_matrix();

        assert_eq!(lorentz_4_4(a, b, &m), 6.0);

        assert_eq!(
            lorentz_4_4(a, b, &m),
            a[0] * b[0] + a[1] * b[1] + a[2] * b[2] - a[3] * b[3]
        );

        assert_eq!(
            lorentz_4_4(a, a, &m),
            a[0].powi(2) + a[1].powi(2) + a[2].powi(2) - a[3].powi(2)
        );
    }
}
