use anise::constants::SPEED_OF_LIGHT_KM_S;

/// Earth angular velocity, in WGS84 frame rad/s
pub const EARTH_ANGULAR_VEL_RAD: f64 = 7.2921151467E-5;

/// Earth gravitational constant (m^3 s-2)
pub const EARTH_GRAVITATION_MU_M3_S2: f64 = 3.986004 * 10.0E14;

// /// Sun gravitational constant (m^3 s-2)
// pub const SUN_GRAVITATION_MU_M3_S2: f64 = 1.327124 * 10.0E20;

// /// Earth's moon gravitational constant (m^3 s-2)
// pub const MOON_GRAVITATION_MU_M3_S2: f64 = 4.902 * 10.0E12;

/// Earth equatorial radis (kilometers)
pub const EARTH_EQUATORIAL_RADIUS_KM: f64 = 6378.1366;

// /// Earth equatorial radis (meters)
// pub const EARTH_EQUATORIAL_RADIUS_M: f64 = EARTH_EQUATORIAL_RADIUS_KM * 1.0E3;

// /// WGS84 Earth Frame Ellipsoid semi-major axis
// pub const EARTH_SEMI_MAJOR_AXIS_WGS84: f64 = 6378137.0_f64;

// /// Love degree^2 term
// pub const LOVE_DEGREE2: f64 = 0.6078;
//
// /// Shida degree^2 term
// pub const SHIDA_DEGREE2: f64 = 0.0847;

/// Speed of light in m.s⁻¹
pub const SPEED_OF_LIGHT_M_S: f64 = SPEED_OF_LIGHT_KM_S * 1000.0;
