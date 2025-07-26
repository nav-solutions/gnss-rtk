use crate::{
    candidate::differences::Differences,
    constants::{
        EARTH_GRAVITATION_MU_M3_S2, // EARTH_SEMI_MAJOR_AXIS_WGS84,
        SPEED_OF_LIGHT_M_S,
    },
    navigation::state::State,
    pool::Pool,
    prelude::{
        // Almanac,
        Candidate,
        EnvironmentalBias,
        EphemerisSource,
        Error,
        OrbitSource,
        SpacebornBias,
        Vector3,
        // EARTH_J2000,
        SUN_J2000,
    },
};

use std::cmp::Ordering;

use log::{debug, error, info};

use hifitime::Unit;

use anise::errors::AlmanacResult;

impl<EPH: EphemerisSource, ORB: OrbitSource, EB: EnvironmentalBias, SB: SpacebornBias>
    Pool<EPH, ORB, EB, SB>
{
    /// Apply Post fit algorithms
    pub fn post_fit(&mut self, name: &str, state: &State) -> AlmanacResult<()> {
        self.post_fit_attitudes(name, state);
        self.post_fit_velocities(name);
        self.post_fit_eclipse(name);
        self.post_fit_biases(name, state);

        Ok(())
    }

    // /// Post fit phase windup correction
    // fn post_fit_phase_windup(&mut self, almanac: &Almanac, _state: &State) -> AlmanacResult<()> {
    //     let epoch = self.inner[0].epoch;

    //     let earth_sun = almanac.transform(SUN_J2000, EARTH_J2000, epoch, None)?;

    //     let r_sun = Vector3::new(
    //         earth_sun.radius_km.x * 1.0E3,
    //         earth_sun.radius_km.y * 1.0E3,
    //         earth_sun.radius_km.z * 1.0E3,
    //     );

    //     // for cd in self.inner.iter_mut() {
    //     // let prev_correction = self
    //     //     .past
    //     //     .iter()
    //     //     .filter_map(|past| {
    //     //         if past.sv == cd.sv {
    //     //             Some(past.windup)
    //     //         } else {
    //     //             None
    //     //         }
    //     //     })
    //     //     .reduce(|k, _| k);

    //     // cd.phase_windup_correction(state, r_sun, prev_correction);
    //     // }

    //     Ok(())
    // }

    /// Apply Attitudes Post fit
    fn post_fit_attitudes(&mut self, name: &str, state: &State) {
        let rx_orbit = state.to_orbit(self.earth_cef);

        self.inner.retain_mut(
            |cd| match cd.orbital_attitude_fixup(&self.almanac, rx_orbit) {
                Ok(_) => {
                    let elevation_deg = cd.elevation_deg.unwrap();
                    if elevation_deg < 0.0 {
                        error!(
                            "{}({}) {} - invalid negative elevation. Invalid input data!",
                            cd.epoch, cd.sv, name,
                        );
                        false
                    } else {
                        true
                    }
                },
                Err(e) => {
                    error!("{}({}) {} orbital fixup: {}", state.epoch, cd.sv, name, e);
                    false
                },
            },
        );

        let min_elev_deg = self.cfg.min_sv_elev.unwrap_or(0.0);
        let min_azim_deg = self.cfg.min_sv_azim.unwrap_or(0.0);
        let max_azim_deg = self.cfg.max_sv_azim.unwrap_or(360.0);

        self.inner.retain(|cd| {
            if let Some((elev, azim)) = cd.attitude() {
                if elev < min_elev_deg {
                    info!("{}({}) - rejected (below elevation mask)", cd.epoch, cd.sv);
                    false
                } else if azim < min_azim_deg {
                    info!("{}({}) - rejected (below azimuth mask)", cd.epoch, cd.sv);
                    false
                } else if azim > max_azim_deg {
                    info!("{}({}) - rejected (above azimuth mask)", cd.epoch, cd.sv);
                    false
                } else {
                    debug!(
                        "{}({}) - elev={:.3}° azim={:.3}°",
                        cd.epoch, cd.sv, elev, azim
                    );
                    true
                }
            } else {
                true
            }
        });
    }

    /// Velocities fit
    fn post_fit_velocities(&mut self, name: &str) {
        // let mu = EARTH_GRAVITATION_MU_M3_S2;
        // let w_e = EARTH_SEMI_MAJOR_AXIS_WGS84;

        for cd in self.inner.iter_mut() {
            let time_of_flight = cd.signal_time_of_flight();

            if let Some(sv_orbit) = &mut cd.orbit {
                if let Some(past_elected) = self.past.iter().find(|elected| elected.sv == cd.sv) {
                    let pos_vel_km_s = sv_orbit.to_cartesian_pos_vel();
                    let past_time_of_flight = past_elected.signal_time_of_flight();

                    let dt_s = (time_of_flight - past_time_of_flight).to_seconds();
                    let past_orbit = past_elected.orbit.unwrap();
                    let pos_vel_z1_km_s = past_orbit.to_cartesian_pos_vel();

                    let vel_km_s = (
                        (pos_vel_km_s[0] - pos_vel_z1_km_s[0]) / dt_s,
                        (pos_vel_km_s[1] - pos_vel_z1_km_s[1]) / dt_s,
                        (pos_vel_km_s[2] - pos_vel_z1_km_s[2]) / dt_s,
                    );

                    debug!(
                        "{}({}) {} inst. velocity {:?} km/s",
                        cd.epoch, cd.sv, name, vel_km_s
                    );

                    *sv_orbit = sv_orbit
                        .with_velocity_km_s(Vector3::new(vel_km_s.0, vel_km_s.1, vel_km_s.2));

                    if self.cfg.modeling.relativistic_clock_bias
                        && cd.clock_corr.needs_relativistic_correction
                    {
                        let pos_m = (
                            pos_vel_km_s[0] * 1.0E3,
                            pos_vel_km_s[1] * 1.0E3,
                            pos_vel_km_s[2] * 1.0E3,
                        );

                        let vel_m_s = (vel_km_s.0 * 1.0E3, vel_km_s.1 * 1.0E3, vel_km_s.2 * 1.0E3);

                        let r_v_sat =
                            pos_m.0 * vel_m_s.0 + pos_m.1 * vel_m_s.1 + pos_m.2 * vel_m_s.2;

                        let bias =
                            -2.0 * r_v_sat / SPEED_OF_LIGHT_M_S / SPEED_OF_LIGHT_M_S * Unit::Second;

                        // let ea_deg = sv_orbit.ea_deg().map_err(Error::Physics)?;

                        // let ea_rad = ea_deg.to_radians();
                        // let gm = (w_e * mu).sqrt();
                        // let ecc = sv_orbit.ecc().map_err(Error::Physics)?;

                        // let bias = -2.0_f64 * ecc * ea_rad.sin() * gm
                        //     / SPEED_OF_LIGHT_M_S
                        //     / SPEED_OF_LIGHT_M_S
                        //     * Unit::Second;

                        debug!(
                            "{}({}) {} : relativistic clock bias: {}",
                            cd.epoch, cd.sv, name, bias
                        );
                        // clock_corr.duration += bias;
                    } //clockbias
                } //velocity
            }
        }
    }

    fn post_fit_eclipse(&mut self, name: &str) {
        self.inner.retain(|cd| {
            let orbit = cd.orbit.unwrap();
            match self
                .almanac
                .occultation(SUN_J2000, self.earth_cef, orbit, None)
            {
                Ok(occultation) => {
                    let retained = occultation.percentage < self.cfg.max_eclipse_rate_percent;

                    if !retained {
                        debug!("{}({}) {} : dropped due to eclipse", cd.epoch, cd.sv, name,);
                    }

                    retained
                },
                Err(e) => {
                    error!("(anise) eclipse error: {e}");
                    false
                },
            }
        });
    }

    // fn post_fit_code_smoothing(&mut self) {
    //     for cd in self.inner.iter_mut() {
    //         for sv_observ in cd.observations.iter_mut() {
    //             if sv_observ.carrier.is_l1() {
    //                 let lambda_1 = sv_observ.carrier.wavelength();
    //                 let n_1 = sv_observ.ambiguity.unwrap_or_default();
    //                 if n_1 != 0.0 {
    //                     if let Some(c_n) = &mut sv_observ.pseudo_range_m {
    //                         if let Some(l_n) = sv_observ.phase_range_m {
    //                             *c_n = self.smoother.smoothing(
    //                                 sv_observ.carrier,
    //                                 cd.sv,
    //                                 *c_n,
    //                                 n_1,
    //                                 lambda_1,
    //                                 l_n,
    //                             );
    //                         }
    //                     }
    //                 }
    //             } else {
    //                 let lambda_j = sv_observ.carrier.wavelength();
    //                 let n_j = sv_observ.ambiguity.unwrap_or_default();
    //                 if n_j != 0.0 {
    //                     if let Some(c_n) = &mut sv_observ.pseudo_range_m {
    //                         if let Some(l_n) = sv_observ.phase_range_m {
    //                             *c_n = self.smoother.smoothing(
    //                                 sv_observ.carrier,
    //                                 cd.sv,
    //                                 *c_n,
    //                                 n_j,
    //                                 lambda_j,
    //                                 l_n,
    //                             );
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    fn post_fit_biases(&mut self, name: &str, state: &State) {
        let rcvr_position_ecef_m = state.to_position_ecef_m();

        for cd in self.inner.iter_mut() {
            let rtm = cd
                .to_bias_runtime(rcvr_position_ecef_m, state.lat_long_alt_deg_deg_km)
                .expect("internal error: post-fit while state is still not resolved!");

            let sat_orbit = cd
                .orbit
                .expect("internal error: post-fit while state is still not resolved!");

            let r_sat_m = sat_orbit.to_cartesian_pos_vel() * 1.0E3;

            if self.cfg.modeling.tropospheric_bias {
                let tropo = self.env_bias.troposphere_bias_m(&rtm);
                debug!("{}({}) {} - tropod={:.3}m", cd.epoch, cd.sv, name, tropo);
                cd.tropod = tropo;
            }

            if self.cfg.modeling.ionospheric_bias {
                let iono = self.env_bias.ionosphere_bias_m(&rtm);
                debug!("{}({}) {} - ionod={:.3}m", cd.epoch, cd.sv, name, iono);
                cd.ionod = iono;
            }

            let pos_ecef_m = state.to_position_ecef_m();
            let (x0_m, y0_m, z0_m) = (pos_ecef_m[0], pos_ecef_m[1], pos_ecef_m[2]);

            let r_0 = (x0_m.powi(2) + y0_m.powi(2) + z0_m.powi(2)).sqrt();
            let r_sat = (r_sat_m[0].powi(2) + r_sat_m[1].powi(2) + r_sat_m[2].powi(2)).sqrt();

            let r_sat_0 = r_0 - r_sat;

            if self.cfg.modeling.relativistic_path_range {
                let dr = 2.0 * EARTH_GRAVITATION_MU_M3_S2 / SPEED_OF_LIGHT_M_S / SPEED_OF_LIGHT_M_S
                    * ((r_sat + r_0 + r_sat_0) / (r_sat + r_0 - r_sat_0)).ln();

                debug!("{}({}) {} - rel. path range={}m", cd.epoch, cd.sv, name, dr);

                cd.relativistic_path_range = dr;
            }

            if self.cfg.modeling.phase_windup {
                // TODO
            }
        }
    }

    /// Runs a postfit SD algorithm, between seld and rhs.
    /// NB: only shared measurements are preserved
    pub fn post_fit_sd(&mut self, pivot: &Candidate) -> Result<(), Error> {
        for cd in self.inner.iter() {
            if cd.sv != pivot.sv {
                self.single_differences
                    .insert(cd.sv, cd.single_difference(pivot));
            }
        }

        // drop pivot
        self.retain(|cd| cd.sv != pivot.sv);

        Ok(())
    }

    /// Select a pivot measurement
    fn post_fit_sd_pivot_election(&self) -> Option<Candidate> {
        self.inner
            .iter()
            .max_by(|meas_a, meas_b| {
                if meas_a.elevation_deg > meas_b.elevation_deg {
                    Ordering::Greater
                } else {
                    Ordering::Less
                }
            })
            .cloned()
    }

    /// Runs the special post-fit prior RTK solving, where self is considered rover
    /// returning DD'ed measurements.
    pub fn rtk_post_fit(&mut self, base: &mut Self) -> Result<Differences, Error> {
        // run SD algorithm on both sites
        let mob_pivot = self.post_fit_sd_pivot_election();

        if mob_pivot.is_none() {
            error!("rover failed to elect pivot satellite");
            return Err(Error::SdPivotSatellite);
        }

        let mob_pivot = mob_pivot.unwrap();

        let base_pivot = base.post_fit_sd_pivot_election();

        if base_pivot.is_none() {
            error!("base failed to elect pivot satellite");
            return Err(Error::SdPivotSatellite);
        }

        let base_pivot = base_pivot.unwrap();

        if base_pivot.sv != mob_pivot.sv {
            error!(
                "{}({}/{}) - pivot sat disagreement - baseline too long!",
                mob_pivot.epoch, mob_pivot.sv, base_pivot.sv
            );
            return Err(Error::RtkBaselineTooLong);
        }

        // SD on both sites
        debug!(
            "{} - using {} as pivot satellite",
            mob_pivot.epoch, mob_pivot.sv
        );

        self.post_fit_sd(&mob_pivot)?;
        base.post_fit_sd(&base_pivot)?;

        let pos_vel = mob_pivot
            .orbit
            .as_ref()
            .expect("internal error: undefined pivot sat state")
            .to_cartesian_pos_vel()
            * 1.0E3;

        self.pivot_position_ecef_m = Some((pos_vel[0], pos_vel[1], pos_vel[2]));

        // DD
        let double_diff = self
            .single_differences
            .double_difference(&base.single_differences);

        // remove pivot from both sites
        self.retain_mut(|cd| cd.sv != mob_pivot.sv);
        base.retain_mut(|cd| cd.sv != mob_pivot.sv);

        for (sat, dd) in double_diff.inner.iter() {
            debug!("{}({}) - DD={}", mob_pivot.epoch, sat, dd);
        }

        Ok(double_diff)
    }
}
