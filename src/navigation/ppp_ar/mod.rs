//! PPP ambiguity resolution (PPP-AR).
//!
//! The float PPP solution (navigation filter, see `navigation::ppp`) provides
//! float ambiguities and their covariance. Fixing them to integers (wide lane
//! from the Melbourne-Wübbena combination, narrow lane with LAMBDA) requires
//! satellite phase biases and is the next milestone.
#[allow(dead_code)]
pub(crate) mod lambda;
