#!/usr/bin/env python3
"""
Download the ANISE data files used by gnss-rtk (de440s.bsp + pck11.pca),
verify their CRC32, and store them where ANISE looks for its cache, or in
any folder of your choice (to be used through GNSS_RTK_ANISE_DIR).

Standard library only (Python >= 3.8).

Examples:
    python3 scripts/fetch_anise_data.py --dest .cache   # repository cache (used by the tests)
    python3 fetch_anise_data.py                       # ANISE cache folder
    python3 fetch_anise_data.py --dest ./anise-data   # local folder
    python3 fetch_anise_data.py --cacert corp-ca.pem  # TLS-intercepting proxy
    python3 fetch_anise_data.py --check               # only verify what is there

Proxy: the usual HTTPS_PROXY / HTTP_PROXY / NO_PROXY environment variables
are honored (urllib).
"""

import argparse
import os
import ssl
import sys
import tempfile
import time
import urllib.request
import zlib
from pathlib import Path

# Name, URL, expected CRC32 (must match the MetaFile entries, ANISE v0.10)
FILES = [
    (
        "de440s.bsp",
        "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440s.bsp",
        0x7286750A,
    ),
    (
        "pck11.pca",
        "http://public-data.nyxspace.com/anise/v0.10/pck11.pca",
        0x1EDB3EAC,
    ),
]

CHUNK = 1 << 20  # 1 MiB


def anise_cache_dir() -> Path:
    """Folder used by ANISE MetaFile (platform_dirs, 'nyx-space/anise', XDG on macOS)."""
    if sys.platform.startswith("win"):
        base = os.environ.get("LOCALAPPDATA") or str(Path.home() / "AppData" / "Local")
        return Path(base) / "nyx-space" / "anise"
    base = os.environ.get("XDG_DATA_HOME") or str(Path.home() / ".local" / "share")
    return Path(base) / "nyx-space" / "anise"


def crc32_of(path: Path) -> int:
    crc = 0
    with open(path, "rb") as fd:
        while True:
            block = fd.read(CHUNK)
            if not block:
                break
            crc = zlib.crc32(block, crc)
    return crc & 0xFFFFFFFF


def human(n: float) -> str:
    for unit in ("B", "KiB", "MiB", "GiB"):
        if n < 1024 or unit == "GiB":
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} GiB"


def download(url: str, dest: Path, expected_crc: int, context, timeout: float) -> None:
    """Streams url into a temporary file next to dest, checks the CRC32, then renames."""
    request = urllib.request.Request(url, headers={"User-Agent": "gnss-rtk-anise-fetch/1.0"})

    fd, tmp_name = tempfile.mkstemp(prefix=dest.name + ".", suffix=".part", dir=dest.parent)
    tmp = Path(tmp_name)

    try:
        with urllib.request.urlopen(request, timeout=timeout, context=context) as resp, os.fdopen(fd, "wb") as out:
            total = int(resp.headers.get("Content-Length") or 0)
            done, crc, t0 = 0, 0, time.time()

            while True:
                block = resp.read(CHUNK)
                if not block:
                    break
                out.write(block)
                crc = zlib.crc32(block, crc)
                done += len(block)

                if total:
                    pct = 100.0 * done / total
                    print(f"\r  {human(done)} / {human(total)} ({pct:5.1f} %)", end="", flush=True)
                else:
                    print(f"\r  {human(done)}", end="", flush=True)

            elapsed = max(time.time() - t0, 1e-3)
            print(f"  [{human(done / elapsed)}/s]")

        crc &= 0xFFFFFFFF
        if crc != expected_crc:
            raise ValueError(f"CRC32 mismatch: got 0x{crc:08x}, expected 0x{expected_crc:08x}")

        os.replace(tmp, dest)  # atomic on the same filesystem
    finally:
        if tmp.exists():
            tmp.unlink()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--dest", type=Path, default=None, help="destination folder (default: ANISE cache folder)")
    parser.add_argument("--cacert", type=Path, default=None, help="extra CA bundle (PEM), e.g. corporate proxy CA")
    parser.add_argument("--timeout", type=float, default=60.0, help="network timeout in seconds (default 60)")
    parser.add_argument("--retries", type=int, default=3, help="attempts per file (default 3)")
    parser.add_argument("--force", action="store_true", help="download even if a valid file is present")
    parser.add_argument("--check", action="store_true", help="only verify the files present, no download")
    args = parser.parse_args()

    dest_dir = (args.dest or anise_cache_dir()).expanduser().resolve()
    dest_dir.mkdir(parents=True, exist_ok=True)
    print(f"Destination: {dest_dir}")

    context = ssl.create_default_context()
    if args.cacert:
        context.load_verify_locations(cafile=str(args.cacert))

    failures = 0

    for name, url, expected in FILES:
        dest = dest_dir / name
        print(f"\n{name}")

        if dest.exists() and not args.force:
            crc = crc32_of(dest)
            if crc == expected:
                print(f"  OK, already present (CRC32 0x{crc:08x})")
                continue
            print(f"  present but CRC32 0x{crc:08x} != 0x{expected:08x}")

        if args.check:
            if not dest.exists():
                print("  missing")
            failures += 1
            continue

        print(f"  downloading {url}")

        for attempt in range(1, args.retries + 1):
            try:
                download(url, dest, expected, context, args.timeout)
                print(f"  OK (CRC32 0x{expected:08x})")
                break
            except Exception as e:  # network, TLS, HTTP or CRC error
                print(f"\n  attempt {attempt}/{args.retries} failed: {e}")
                if isinstance(e, ssl.SSLError) or "CERTIFICATE_VERIFY_FAILED" in str(e):
                    print("  hint: TLS interception? retry with --cacert <your proxy CA.pem>")
                if attempt < args.retries:
                    time.sleep(2 * attempt)
        else:
            failures += 1

    print()
    if failures:
        print(f"{failures} file(s) missing or invalid.")
        return 1

    print("All files present and valid.")
    if args.dest:
        print(f"Use them with:  export GNSS_RTK_ANISE_DIR={dest_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
