//! Builds an ANISE Almanac for complete test suite

use anise::{
    almanac::metaload::{MetaAlmanac, MetaFile},
    constants::frames::EARTH_J2000,
    prelude::{Almanac, Epoch, Frame, Orbit},
};

use std::{
    path::{Path, PathBuf},
    str::FromStr,
};

use crate::navigation::apriori::Apriori;

use log::LevelFilter;
use std::sync::{Once, OnceLock};

static ALMANAC: OnceLock<Almanac> = OnceLock::new();

pub const ANISE_DATA_FILES: [(&str, &str, u32); 2] = [
    (
        "de440s.bsp",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
        0x7286750a,
    ),
    (
        "pck11.pca",
        "http://public-data.nyxspace.com/anise/v0.10/pck11.pca",
        0x1edb3eac,
    ),
];

/// Local folder holding the ANISE data files:
/// `$GNSS_RTK_ANISE_DIR` if defined, `<repo>/.cache` otherwise.
pub fn anise_cache_dir() -> PathBuf {
    match std::env::var("GNSS_RTK_ANISE_DIR") {
        Ok(dir) => PathBuf::from(dir),
        Err(_) => Path::new(env!("CARGO_MANIFEST_DIR")).join(".cache"),
    }
}

fn load_almanac() -> Almanac {
    let dir = anise_cache_dir();

    // Local files are used as is: no network access at all.
    // Missing files are downloaded by ANISE (cached and verified by CRC32
    // in its own data folder).
    let files = ANISE_DATA_FILES
        .iter()
        .map(|(name, remote, crc32)| {
            let local = dir.join(name);
            if local.is_file() {
                MetaFile {
                    uri: local.to_string_lossy().to_string(),
                    crc32: None,
                }
            } else {
                MetaFile {
                    uri: remote.to_string(),
                    crc32: Some(*crc32),
                }
            }
        })
        .collect();

    MetaAlmanac { files }.process(true).unwrap_or_else(|e| {
        panic!(
            "Failed to build test Almanac: {e}\n\
            Place de440s.bsp and pck11.pca in {} 
            for example with: python3 scripts/fetch_anise_data.py --dest {}",
            dir.display(),
            dir.display(),
        )
    })
}

pub fn almanac() -> Almanac {
    ALMANAC.get_or_init(load_almanac).clone()
}

pub fn earth_frame() -> Frame {
    let alm = almanac();
    alm
        // .unwrap_or_else(|e| panic!("Failed to build test Almanac: {e}"))
        .frame_from_uid(EARTH_J2000)
        .unwrap_or_else(|e| panic!("Failed to build test EARTH-J2000 frame: {e}"))
}
