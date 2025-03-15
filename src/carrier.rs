use crate::prelude::Error;
use nyx::cosmic::SPEED_OF_LIGHT_M_S;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Copy, Default, PartialEq, PartialOrd, Eq, Ord, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Carrier {
    /// L1 (GPS/QZSS/SBAS) same frequency as E1 and B1aB1c
    #[default]
    L1,
    /// L2 (GPS/QZSS)
    L2,
    /// L5 (GPS/QZSS/SBAS) same frequency as E5A and B2A
    L5,
    /// L6 (GPS/QZSS) same frequency as E6
    L6,
    /// E1 (Galileo)
    E1,
    /// E5 (Galileo) same frequency as B2
    E5,
    /// E5A (Galileo) same frequency as L5
    E5A,
    /// E5B (Galileo) same frequency as B2iB2b
    E5B,
    /// E6 (Galileo) same frequency as L6
    E6,
    /// B1aB1c (BDS) same frequency as L1
    B1aB1c,
    /// B1I (BDS)
    B1I,
    /// B2I/B2B (BDS) same frequency as E5b
    B2iB2b,
    /// B2 (BDS) same frequency as E5
    B2,
    /// B2A (BDS) same frequency as L5 and E5A
    B2A,
    /// B3 (BDS)
    B3,
}

impl std::fmt::Display for Carrier {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        match self {
            Self::L1 => write!(f, "L1"),
            Self::L2 => write!(f, "L2"),
            Self::L5 => write!(f, "L5"),
            Self::L6 => write!(f, "L6"),
            Self::E1 => write!(f, "E1"),
            Self::E5 => write!(f, "E5"),
            Self::E5A => write!(f, "E5A"),
            Self::E5B => write!(f, "E5B"),
            Self::E6 => write!(f, "E6"),
            Self::B1I => write!(f, "B1I"),
            Self::B1aB1c => write!(f, "B1A/B1C"),
            Self::B2iB2b => write!(f, "B2I/B2B"),
            Self::B2 => write!(f, "B2"),
            Self::B3 => write!(f, "B3"),
            Self::B2A => write!(f, "B2A"),
        }
    }
}

impl std::str::FromStr for Carrier {
    type Err = Error;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let trimmed = s.trim();
        if trimmed.eq("L1") {
            Ok(Self::L1)
        } else if trimmed.eq("L2") {
            Ok(Self::L2)
        } else if trimmed.eq("L5") {
            Ok(Self::L5)
        } else if trimmed.eq("L6") {
            Ok(Self::L6)
        } else if trimmed.eq("E1") {
            Ok(Self::E1)
        } else if trimmed.eq("E5") {
            Ok(Self::E5)
        } else if trimmed.eq("E6") {
            Ok(Self::E6)
        } else if trimmed.eq("E5A") {
            Ok(Self::E5A)
        } else if trimmed.eq("E5B") {
            Ok(Self::E5B)
        } else if trimmed.eq("B1I") {
            Ok(Self::B1I)
        } else if trimmed.eq("B2") {
            Ok(Self::B2)
        } else if trimmed.eq("B3") {
            Ok(Self::B3)
        } else if trimmed.eq("B2A") {
            Ok(Self::B2A)
        } else if trimmed.contains("B1A") {
            Ok(Self::B1aB1c)
        } else if trimmed.contains("B1C") {
            Ok(Self::B1aB1c)
        } else if trimmed.contains("B2I") {
            Ok(Self::B2iB2b)
        } else if trimmed.contains("B2B") {
            Ok(Self::B2iB2b)
        } else {
            Err(Error::InvalidFrequency)
        }
    }
}

impl Carrier {
    pub fn frequency(&self) -> f64 {
        match self {
            Self::L1 | Self::E1 | Self::B1aB1c => 1575.42E6_f64,
            Self::L2 => 1227.60E6_f64,
            Self::L5 | Self::E5A | Self::B2A => 1176.45E6_f64,
            Self::E5 | Self::B2 => 1191.795E6_f64,
            Self::L6 | Self::E6 => 1278.750E6_f64,
            Self::B3 => 1268.52E6_f64,
            Self::E5B | Self::B2iB2b => 1207.14E6_f64,
            Self::B1I => 1561.098E6_f64,
        }
    }
    pub fn wavelength(&self) -> f64 {
        SPEED_OF_LIGHT_M_S / self.frequency()
    }
}

/// Signal used in [PVTSolution] resolution
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Signal {
    Single(Carrier),
    Dual((Carrier, Carrier)),
}

impl Default for Signal {
    fn default() -> Self {
        Self::Single(Default::default())
    }
}

impl std::fmt::Display for Signal {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        match self {
            Self::Single(carrier) => write!(f, "{}", carrier),
            Self::Dual((lhs, rhs)) => write!(f, "{}/{}", rhs, lhs),
        }
    }
}

impl Signal {
    pub(crate) fn single(carrier: Carrier) -> Self {
        Self::Single(carrier)
    }
    pub(crate) fn dual(lhs: Carrier, rhs: Carrier) -> Self {
        Self::Dual((lhs, rhs))
    }
}
