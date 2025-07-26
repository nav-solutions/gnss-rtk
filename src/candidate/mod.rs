//! Position solving candidate
use hifitime::Unit;
use log::debug;

// use itertools::Itertools;

use crate::{
    bias::spaceborn::SatelliteClockCorrection,
    constants::SPEED_OF_LIGHT_M_S,
    navigation::state::State,
    prelude::{Almanac, Config, Duration, Epoch, Error, Orbit, Vector3, SV},
};

use anise::errors::AlmanacResult;

mod bias;
mod ppp;
mod rtk;
mod signal;

pub(crate) mod combination;
pub(crate) mod differences;
pub(crate) mod single_diff;

#[cfg(test)]
mod tests;

pub use crate::{candidate::signal::Observation, prelude::Carrier};

/// Position solving candidate
#[derive(Clone, Debug)]
pub struct Candidate {
    /// [SV]
    pub sv: SV,

    /// Sampling [Epoch]
    pub epoch: Epoch,

    /// Transmission [Epoch]
    pub(crate) tx_epoch: Epoch,

    /// [Orbit]al state
    pub(crate) orbit: Option<Orbit>,

    /// elevation at reception time
    pub(crate) elevation_deg: Option<f64>,

    /// azimuth at reception time
    pub(crate) azimuth_deg: Option<f64>,

    /// Total group delay.
    pub(crate) tgd: Duration,

    /// Troposphere delay (meters)
    pub(crate) tropod: f64,

    /// Ionosphere delay (meters)
    pub(crate) ionod: f64,

    // /// Phase wind up in cycles
    // pub(crate) windup: f64,
    /// [SatelliteClockCorrection]
    pub(crate) clock_corr: SatelliteClockCorrection,

    /// Signal [Observation]s
    pub(crate) observations: Vec<Observation>,

    /// Possible time system correction
    pub(crate) system_correction: Option<Duration>,

    /// Estimated relativistic path range
    pub(crate) relativistic_path_range: f64,
}

impl Candidate {
    /// Basic candidate definition, to propose to the navigation solver. Each candidate is to be
    /// This is the most simplistic definition (bare minimum).
    /// It is certainly not enough for PPP navigation and will require
    /// you provide more information (see other customization methods),
    /// especially if you want to achieve accurate results.
    /// ## Input
    /// - sv: [SV] Identity
    /// - epoch: sampling [Epoch]
    /// - observations: provide signals observations.
    ///   You have to provide observations that match your navigation method.
    pub fn new(sv: SV, epoch: Epoch, observations: Vec<Observation>) -> Self {
        Self {
            sv,
            epoch,
            ionod: 0.0,
            tropod: 0.0,
            observations,
            tx_epoch: epoch,
            tgd: Duration::ZERO,
            orbit: Default::default(),
            // windup: Default::default(),
            azimuth_deg: Default::default(),
            clock_corr: Default::default(),
            elevation_deg: Default::default(),
            system_correction: Default::default(),
            relativistic_path_range: Default::default(),
        }
    }

    // /// Designs a measured frequency iterator
    // fn frequencies_iter(&self) -> Box<dyn Iterator<Item = Carrier> + '_> {
    //     Box::new(self.observations.iter().map(|obs| obs.carrier).unique())
    // }

    /// Update pseudo range observation (in meters) for this frequency.
    pub fn set_pseudo_range_m(&mut self, carrier: Carrier, pr_m: f64) {
        if let Some(observation) = self
            .observations
            .iter_mut()
            .filter_map(|obs| {
                if obs.carrier == carrier && obs.pseudo_range_m.is_some() {
                    Some(obs)
                } else {
                    None
                }
            })
            .reduce(|k, _| k)
        {
            observation.pseudo_range_m = Some(pr_m);
        } else {
            self.observations
                .push(Observation::pseudo_range(carrier, pr_m, None));
        }
    }

    /// Update with ambiguous range observation (in meters) for this frequency.
    pub fn set_ambiguous_phase_range_m(&mut self, carrier: Carrier, pr_m: f64) {
        if let Some(observation) = self
            .observations
            .iter_mut()
            .filter_map(|obs| {
                if obs.carrier == carrier && obs.phase_range_m.is_some() {
                    Some(obs)
                } else {
                    None
                }
            })
            .reduce(|k, _| k)
        {
            observation.phase_range_m = Some(pr_m);
        } else {
            self.observations
                .push(Observation::ambiguous_phase_range(carrier, pr_m, None));
        }
    }

    /// Computes phase windup correction term.
    pub(crate) fn phase_windup_correction(
        &mut self,
        _rx_state: &State,
        _r_sun: Vector3<f64>,
        _past_correction: Option<f64>,
    ) {
        // let sv_state = self.orbit.unwrap_or_else(|| {
        //     panic!("internal error: phase windup while state is not fully resolved");
        // });

        // let r_sv = sv_state.to_cartesian_pos_vel() * 1.0E3;

        // // todo self.yaw_attitude();

        // let r_rx = rx_state.to_position_ecef_m();
        // let r_sv = Vector3::new(r_sv[0], r_sv[1], r_sv[2]);
        // let r_rr_rs = r_rx - r_sv;
        // let e_r_rs = r_rr_rs.norm();

        // self.windup = 0.0; // TODO
    }

    /// Computes signal transmission instant, as [Epoch].
    pub(crate) fn transmission_time(&mut self, name: &str, cfg: &Config) -> Result<(), Error> {
        let mut t_tx = self.epoch;

        let (_, pr) = self.best_snr_range_m().ok_or(Error::MissingPseudoRange)?;

        t_tx -= pr / SPEED_OF_LIGHT_M_S * Unit::Second;

        if cfg.modeling.sv_total_group_delay {
            t_tx -= self.tgd;
        }

        if cfg.modeling.sv_clock_bias {
            t_tx -= self.clock_corr.duration;
        }

        assert!(
            t_tx < self.epoch,
            "Physical non sense - {} rx={} prior tx={}",
            name,
            self.epoch,
            t_tx
        );

        self.tx_epoch = t_tx;

        debug!(
            "{}({}) {} - time of flight: {}",
            self.sv,
            self.epoch,
            name,
            self.signal_time_of_flight()
        );

        Ok(())
    }

    pub(crate) fn signal_time_of_flight(&self) -> Duration {
        self.epoch - self.tx_epoch
    }

    /// Fix [Orbit]al attitude.
    pub(crate) fn orbital_attitude_fixup(
        &mut self,
        almanac: &Almanac,
        rx_orbit: Orbit,
    ) -> AlmanacResult<()> {
        let orbit = self
            .orbit
            .expect("internal error: undefined orbital state (badop)");

        let elazrg = almanac.azimuth_elevation_range_sez(orbit, rx_orbit, None, None)?;

        self.azimuth_deg = Some(elazrg.azimuth_deg);
        self.elevation_deg = Some(elazrg.elevation_deg);
        Ok(())
    }

    /// Returns (elevation, azimuth) in decimal degrees.
    pub(crate) fn attitude(&self) -> Option<(f64, f64)> {
        let el = self.elevation_deg?;
        let az = self.azimuth_deg?;
        Some((el, az))
    }

    // /// Creates [Candidate] with elevation degree
    // pub(crate) fn with_elevation_deg(&self, el: f64) -> Self {
    //     let mut s = self.clone();
    //     s.elevation_deg = Some(el);
    //     s
    // }

    // /// Creates [Candidate] with azimuth degree
    // pub(crate) fn with_azimuth_deg(&self, az: f64) -> Self {
    //     let mut s = self.clone();
    //     s.azimuth_deg = Some(az);
    //     s
    // }

    #[cfg(test)]
    pub fn set_orbit(&mut self, orbit: Orbit) {
        self.orbit = Some(orbit);
    }
}

#[cfg(test)]
mod test {
    use crate::prelude::{Candidate, Carrier, Epoch, Observation, SV};

    #[test]
    fn cpp_compatibility() {
        {
            let (observations, cpp_compatible) = (
                vec![
                    Observation {
                        snr_dbhz: Some(1.0),
                        pseudo_range_m: Some(1.0),
                        phase_range_m: Some(2.0),
                        ambiguity: None,
                        doppler: None,
                        carrier: Carrier::L1,
                    },
                    Observation {
                        snr_dbhz: Some(1.0),
                        pseudo_range_m: Some(2.0),
                        phase_range_m: Some(2.0),
                        ambiguity: None,
                        doppler: None,
                        carrier: Carrier::L5,
                    },
                ],
                true,
            );
            let cd = Candidate::new(SV::default(), Epoch::default(), observations);
            assert_eq!(cd.cpp_compatible(), cpp_compatible);
        }
    }
}
