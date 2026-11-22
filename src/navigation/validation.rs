use log::debug;
use nalgebra::{DVector, Vector3};
use nyx::cosmic::SPEED_OF_LIGHT_M_S;
use thiserror::Error;

use crate::{
    navigation::{Input, Output, PVTSolutionType},
    prelude::{Candidate, Config},
};

#[derive(Clone, Debug, PartialEq, Error)]
/// Reason why this solution has been invalidated
pub enum InvalidationCause {
    /// First solution is always invalidated
    #[error("first solution")]
    FirstSolution,
    #[error("gdop {0}: limit exceeded")]
    GDOPOutlier(f64),
    #[error("tdop limit exceeded {0}")]
    TDOPOutlier(f64),
    #[error("innovation outlier |{0}|")]
    InnovationOutlier(f64),
    #[error("coderes limit exceeded {0}")]
    CodeResidual(f64),
}

pub(crate) struct Validator {
    gdop: f64,
    tdop: f64,
    residuals: DVector<f64>,
}

impl Validator {
    pub fn new(apriori: Vector3<f64>, pool: &[Candidate], input: &Input, output: &Output) -> Self {
        let gdop = output.gdop;
        let tdop = output.tdop;
        let mut residuals = DVector::<f64>::zeros(pool.len());

        for (idx, cd) in pool.iter().enumerate() {
            let sv = input
                .sv
                .iter()
                .filter_map(|(sv, data)| if *sv == cd.sv { Some(data) } else { None })
                .reduce(|k, _| k)
                .unwrap();

            let pr = cd.prefered_pseudorange().unwrap().pseudo.unwrap();

            let x = output.state.estimate();

            let (x, y, z, dt) = (
                apriori[0] + x[0],
                apriori[1] + x[1],
                apriori[2] + x[2],
                x[3] / SPEED_OF_LIGHT_M_S,
            );

            let sv_orbit = cd.orbit.unwrap();
            let (sv_x, sv_y, sv_z) = (
                sv_orbit.radius_km.x * 1.0E3,
                sv_orbit.radius_km.y * 1.0E3,
                sv_orbit.radius_km.z * 1.0E3,
            );

            let rho = ((sv_x - x).powi(2) + (sv_y - y).powi(2) + (sv_z - z).powi(2)).sqrt();

            let dt = cd.clock_corr.unwrap().seconds() - dt;

            residuals[idx] = pr;
            residuals[idx] -= rho;
            residuals[idx] += dt * SPEED_OF_LIGHT_M_S;
            residuals[idx] -= sv.tropo_bias.unwrap_or_default();
            residuals[idx] -= sv.iono_bias.unwrap_or_default().value();
            residuals[idx] /= input.w[(idx, idx)];
            debug!(
                "{} ({}): coderes={}/w={}",
                cd.t,
                cd.sv,
                residuals[idx],
                input.w[(idx, idx)]
            );
        }
        Self {
            residuals,
            gdop,
            tdop,
        }
    }
    /*
     * Solution validation process
     */
    pub fn validate(&self, cfg: &Config) -> Result<(), InvalidationCause> {
        if cfg.sol_type != PVTSolutionType::TimeOnly {
            // Other geometry criteria apply
            if let Some(max_gdop) = cfg.solver.gdop_threshold {
                if self.gdop > max_gdop {
                    return Err(InvalidationCause::GDOPOutlier(self.gdop));
                }
            }
            if let Some(max_tdop) = cfg.solver.tdop_threshold {
                if self.tdop > max_tdop {
                    return Err(InvalidationCause::TDOPOutlier(self.tdop));
                }
            }
        }
        Ok(())
    }
}
