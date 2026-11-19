#!/usr/bin/env python3
"""
Extracts frozen test data for gnss-rtk from RINEX / SP3 / CLK files.

Inputs (nav-solutions/data, 2020-06-25):
  - observation RINEX V3 (decompress CRINEX first, e.g. with crx2rnx)
  - SP3 orbits (GPS, Galileo), the day before is useful for interpolation
  - CLK clocks (30s)
  - NAV RINEX V3 (BeiDou broadcast ephemerides: no BDS in the SP3)

Outputs (text, one record per line, see the headers):
  - <station>.obs : raw observations (code m, phase m, SNR dB-Hz) per epoch / satellite
  - orbits.dat    : orbital samples (ECEF km) to be interpolated at any instant
  - clocks.dat    : satellite clocks at the observation epochs (s), relativistic
                    correction (s) and group delay to apply (s)

Standard library only. Example:

  crx2rnx ESBC00DNK_R_20201770000_01D_30S_MO.crx.gz -o ESBC.rnx
  crx2rnx MOJN00DNK_R_20201770000_01D_30S_MO.crx.gz -o MOJN.rnx
  python3 scripts/extract_frozen_data.py \\
      --obs ESBC.rnx --obs MOJN.rnx \\
      --sp3 GRG0MGXFIN_20201760000_01D_15M_ORB.SP3.gz \\
      --sp3 GRG0MGXFIN_20201770000_01D_15M_ORB.SP3.gz \\
      --clk GRG0MGXFIN_20201770000_01D_30S_CLK.CLK.gz \\
      --nav ESBC00DNK_R_20201770000_01D_MN.rnx.gz \\
      --start "2020-06-25 00:00:00" --duration 3600 --interval 30 \\
      --out src/tests/frozen/data/dense
"""

import argparse
import gzip
import math
import os
from datetime import datetime, timedelta

C = 299_792_458.0
GPS_EPOCH = datetime(1980, 1, 6)
BDT_MINUS_GPST_S = -14.0

# BeiDou ICD constants
BDS_GM = 3.986004418e14
BDS_OMEGA_E = 7.2921150e-5
REL_F = -4.442807309e-10

FREQ_HZ = {
    "L1": 1575.42e6,
    "L2": 1227.60e6,
    "L5": 1176.45e6,
    "B1": 1561.098e6,
    "B3": 1268.52e6,
}

# Signals per constellation: (solver carrier, code candidates, phase candidates, snr candidates).
# GPS & Galileo signals match the GRG clock datum (iono free P1/P2 and E1/E5a).
SIGNALS = {
    "G": [
        ("L1", ["C1W", "C1C"], ["L1C", "L1W"], ["S1W", "S1C"]),
        ("L2", ["C2W", "C2L"], ["L2W", "L2L"], ["S2W", "S2L"]),
    ],
    "E": [
        ("L1", ["C1C"], ["L1C"], ["S1C"]),
        ("L5", ["C5Q"], ["L5Q"], ["S5Q"]),
    ],
    "C": [
        ("B1", ["C2I"], ["L2I"], ["S2I"]),
        ("B3", ["C6I"], ["L6I"], ["S6I"]),
    ],
}


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")


def fnum(s):
    s = s.strip().replace("D", "E").replace("d", "e")
    return float(s) if s else None


def parse_time(s):
    return datetime.strptime(s, "%Y-%m-%d %H:%M:%S")


# ---------------------------------------------------------------- observations
def read_obs(path, t0, t1, interval_s, systems):
    header = {"types": {}}
    with opener(path) as fd:
        cur_sys = None
        for line in fd:
            label = line[60:].strip()
            if label == "END OF HEADER":
                break
            if label == "MARKER NAME":
                header["marker"] = line[:60].strip()
            elif label == "APPROX POSITION XYZ":
                header["xyz"] = [float(v) for v in line[:42].split()]
            elif label == "ANTENNA: DELTA H/E/N":
                header["hen"] = [float(v) for v in line[:42].split()]
            elif label == "SYS / # / OBS TYPES":
                if line[0] != " ":
                    cur_sys = line[0]
                    header["types"][cur_sys] = []
                header["types"][cur_sys] += line[7:60].split()

        records = []
        keep = False
        t = None
        for line in fd:
            if line.startswith(">"):
                y, mo, d, h, mi = (int(v) for v in line[2:18].split())
                sec = float(line[18:29])
                t = datetime(y, mo, d, h, mi) + timedelta(seconds=sec)
                flag = int(line[29:32])
                dt = (t - t0).total_seconds()
                keep = flag == 0 and t0 <= t <= t1 and abs(dt / interval_s - round(dt / interval_s)) < 1e-6
                continue
            if not keep or line[0] not in systems:
                continue
            sys_ = line[0]
            sv = line[:3].replace(" ", "0")
            types = header["types"][sys_]
            values = {}
            for i, code in enumerate(types):
                field = line[3 + 16 * i: 3 + 16 * i + 14]
                v = fnum(field) if field.strip() else None
                if v is not None:
                    values[code] = v
            sigs = []
            for carrier, codes, phases, snrs in SIGNALS[sys_]:
                c_m = next((values[k] for k in codes if k in values), None)
                if c_m is None:
                    continue
                l_cyc = next((values[k] for k in phases if k in values), None)
                s_db = next((values[k] for k in snrs if k in values), None)
                l_m = l_cyc * C / FREQ_HZ[carrier] if l_cyc is not None else None
                sigs.append((carrier, c_m, l_m, s_db))
            # dual frequency only (all methods can then be tested on the same data)
            if len(sigs) == 2:
                records.append(((t - t0).total_seconds(), sv, sigs))
    return header, records


# ---------------------------------------------------------------- SP3 / CLK
def read_sp3(paths, systems):
    samples = {}
    for path in paths:
        with opener(path) as fd:
            t = None
            for line in fd:
                if line.startswith("*  "):
                    y, mo, d, h, mi = (int(v) for v in line[3:20].split())
                    t = datetime(y, mo, d, h, mi) + timedelta(seconds=float(line[20:31]))
                elif line.startswith("P") and line[1] in systems and t is not None:
                    sv = line[1:4].replace(" ", "0")
                    x, y_, z = (float(line[4 + 14 * i: 18 + 14 * i]) for i in range(3))
                    if abs(x) < 1e-6 and abs(y_) < 1e-6:
                        continue
                    samples.setdefault(sv, {})[t] = (x, y_, z)
    return samples


def read_clk(path, systems, wanted_epochs):
    clocks = {}
    with opener(path) as fd:
        for line in fd:
            if not line.startswith("AS ") or line[3] not in systems:
                continue
            sv = line[3:6].replace(" ", "0")
            parts = line[8:].split()
            y, mo, d, h, mi = (int(v) for v in parts[:5])
            t = datetime(y, mo, d, h, mi) + timedelta(seconds=float(parts[5]))
            if t in wanted_epochs:
                clocks[(sv, t)] = float(parts[7])
    return clocks


def lagrange(samples, t, order=10):
    """Lagrange interpolation of {datetime: (x,y,z)} at t using order+1 nearest samples."""
    keys = sorted(samples, key=lambda k: abs((k - t).total_seconds()))[: order + 1]
    if len(keys) < order + 1:
        return None
    ts = [(k - t).total_seconds() for k in keys]
    out = [0.0, 0.0, 0.0]
    for i, k in enumerate(keys):
        li = 1.0
        for j in range(len(keys)):
            if j != i:
                li *= (0.0 - ts[j]) / (ts[i] - ts[j])
        for a in range(3):
            out[a] += li * samples[k][a]
    return out


def relativistic_s(samples, t):
    """-2 r.v / c² from SP3 samples (km), at t."""
    p = lagrange(samples, t)
    pm = lagrange(samples, t - timedelta(seconds=0.5))
    pp = lagrange(samples, t + timedelta(seconds=0.5))
    if p is None or pm is None or pp is None:
        return None
    r = [v * 1e3 for v in p]
    v = [(pp[a] - pm[a]) * 1e3 for a in range(3)]
    return -2.0 * sum(r[a] * v[a] for a in range(3)) / C**2


# ---------------------------------------------------------------- BeiDou broadcast
def read_bds_nav(path):
    ephs = {}
    with opener(path) as fd:
        for line in fd:
            if line[60:].strip() == "END OF HEADER":
                break
        lines = fd.readlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line or line[0] == " ":
            i += 1
            continue
        sys_ = line[0]
        n = {"G": 8, "E": 8, "C": 8, "J": 8, "I": 8, "R": 4, "S": 4}.get(sys_, 8)
        block = lines[i: i + n]
        i += n
        if sys_ != "C":
            continue
        sv = line[:3].replace(" ", "0")
        y, mo, d, h, mi, s = (int(v) for v in line[4:23].split())
        toc = datetime(y, mo, d, h, mi, s)  # BDT

        def f(row, col):
            if row == 0:
                return fnum(block[0][23 + 19 * (col - 1): 23 + 19 * col])
            return fnum(block[row][4 + 19 * col: 4 + 19 * (col + 1)])

        e = {
            "toc": toc,
            "af0": f(0, 1), "af1": f(0, 2), "af2": f(0, 3),
            "crs": f(1, 1), "dn": f(1, 2), "m0": f(1, 3),
            "cuc": f(2, 0), "ecc": f(2, 1), "cus": f(2, 2), "sqrta": f(2, 3),
            "toe": f(3, 0), "cic": f(3, 1), "omg0": f(3, 2), "cis": f(3, 3),
            "i0": f(4, 0), "crc": f(4, 1), "omg": f(4, 2), "omgd": f(4, 3),
            "idot": f(5, 0), "week": f(5, 2),
            "tgd1": f(6, 2),
        }
        ephs.setdefault(sv, []).append(e)
    return ephs


def bds_sow(t_gpst):
    t_bdt = t_gpst + timedelta(seconds=BDT_MINUS_GPST_S)
    gps_s = (t_bdt - GPS_EPOCH).total_seconds() - 1356 * 604800
    return t_bdt, gps_s % 604800


def bds_select(ephs, t_gpst):
    t_bdt, _ = bds_sow(t_gpst)
    return min(ephs, key=lambda e: abs((e["toc"] - t_bdt).total_seconds()))


def bds_state(e, sv, t_gpst):
    """ECEF position (km), clock (s), relativistic correction (s), BDS ICD."""
    t_bdt, sow = bds_sow(t_gpst)
    a = e["sqrta"] ** 2
    tk = sow - e["toe"]
    if tk > 302400:
        tk -= 604800
    elif tk < -302400:
        tk += 604800
    n = math.sqrt(BDS_GM / a**3) + e["dn"]
    m = e["m0"] + n * tk
    ecc = e["ecc"]
    ek = m
    for _ in range(30):
        ek = m + ecc * math.sin(ek)
    v = math.atan2(math.sqrt(1 - ecc**2) * math.sin(ek), math.cos(ek) - ecc)
    phi = v + e["omg"]
    s2, c2 = math.sin(2 * phi), math.cos(2 * phi)
    u = phi + e["cus"] * s2 + e["cuc"] * c2
    r = a * (1 - ecc * math.cos(ek)) + e["crs"] * s2 + e["crc"] * c2
    i = e["i0"] + e["idot"] * tk + e["cis"] * s2 + e["cic"] * c2
    xp, yp = r * math.cos(u), r * math.sin(u)
    prn = int(sv[1:])
    if prn <= 5 or prn >= 59:
        # GEO
        om = e["omg0"] + e["omgd"] * tk - BDS_OMEGA_E * e["toe"]
        xg = xp * math.cos(om) - yp * math.cos(i) * math.sin(om)
        yg = xp * math.sin(om) + yp * math.cos(i) * math.cos(om)
        zg = yp * math.sin(i)
        f5 = math.radians(-5.0)
        x1, y1, z1 = xg, yg * math.cos(f5) + zg * math.sin(f5), -yg * math.sin(f5) + zg * math.cos(f5)
        z_ = BDS_OMEGA_E * tk
        x = x1 * math.cos(z_) + y1 * math.sin(z_)
        y = -x1 * math.sin(z_) + y1 * math.cos(z_)
        z = z1
    else:
        om = e["omg0"] + (e["omgd"] - BDS_OMEGA_E) * tk - BDS_OMEGA_E * e["toe"]
        x = xp * math.cos(om) - yp * math.cos(i) * math.sin(om)
        y = xp * math.sin(om) + yp * math.cos(i) * math.cos(om)
        z = yp * math.sin(i)
    dtc = (t_bdt - e["toc"]).total_seconds()
    clk = e["af0"] + e["af1"] * dtc + e["af2"] * dtc**2
    rel = REL_F * ecc * e["sqrta"] * math.sin(ek)
    return (x / 1e3, y / 1e3, z / 1e3), clk, rel


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--obs", action="append", required=True)
    ap.add_argument("--sp3", action="append", required=True)
    ap.add_argument("--clk", required=True)
    ap.add_argument("--nav", required=True)
    ap.add_argument("--start", required=True, help="GPST, 'YYYY-MM-DD hh:mm:ss'")
    ap.add_argument("--duration", type=float, required=True, help="seconds")
    ap.add_argument("--interval", type=float, required=True, help="seconds")
    ap.add_argument("--systems", default="GEC")
    ap.add_argument("--orbit-margin", type=float, default=3 * 3600.0, help="orbit samples margin (s)")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    t0 = parse_time(args.start)
    t1 = t0 + timedelta(seconds=args.duration)
    systems = set(args.systems)
    os.makedirs(args.out, exist_ok=True)

    stations = [read_obs(path, t0, t1, args.interval, systems) for path in args.obs]

    sp3 = read_sp3(args.sp3, systems - {"C"})
    ephs = read_bds_nav(args.nav) if "C" in systems else {}

    all_epochs = sorted({t0 + timedelta(seconds=dt) for _, records in stations for dt, _, _ in records})
    clk = read_clk(args.clk, systems - {"C"}, set(all_epochs))

    # satellite clock (s), relativistic correction (s), group delay (s) per (sv, epoch);
    # observations without precise products are discarded (like any PPP software).
    products = {}
    for _, records in stations:
        for dt, sv, _ in records:
            t = t0 + timedelta(seconds=dt)
            if (sv, t) in products:
                continue
            if sv[0] == "C":
                if sv not in ephs:
                    continue
                e = bds_select(ephs[sv], t)
                if abs((e["toc"] - (t + timedelta(seconds=BDT_MINUS_GPST_S))).total_seconds()) > 4 * 3600:
                    continue
                _, c_s, rel = bds_state(e, sv, t)
                products[(sv, t)] = (c_s, rel, -e["tgd1"])
            elif (sv, t) in clk and sv in sp3:
                rel = relativistic_s(sp3[sv], t)
                if rel is not None:
                    products[(sv, t)] = (clk[(sv, t)], rel, 0.0)

    epochs = set()
    used_sv = set()

    for header, records in stations:
        name = header["marker"].lower()
        out = os.path.join(args.out, f"{name}.obs")
        n = 0
        with open(out, "w") as fd:
            fd.write(f"# station {header['marker']}\n")
            fd.write("# marker_ecef_m {:.4f} {:.4f} {:.4f}\n".format(*header["xyz"]))
            fd.write("# antenna_hen_m {:.4f} {:.4f} {:.4f}\n".format(*header["hen"]))
            fd.write(f"# t0 {t0:%Y-%m-%dT%H:%M:%S} GPST\n")
            fd.write("# columns: dt_s sv [carrier code_m phase_m snr_dbhz] x2 ('-' when missing)\n")
            for dt, sv, sigs in records:
                t = t0 + timedelta(seconds=dt)
                if (sv, t) not in products:
                    continue
                cols = [f"{dt:g}", sv]
                for carrier, c_m, l_m, s_db in sigs:
                    cols += [
                        carrier,
                        f"{c_m:.3f}",
                        f"{l_m:.3f}" if l_m is not None else "-",
                        f"{s_db:.0f}" if s_db is not None else "-",
                    ]
                fd.write(" ".join(cols) + "\n")
                epochs.add(t)
                used_sv.add(sv)
                n += 1
        print(f"{out}: {n} records ({len(records) - n} discarded: no precise products)")

    epochs = sorted(epochs)

    # orbits: SP3 (GPS, Galileo) samples, BDS broadcast samples every 15'
    o0 = t0 - timedelta(seconds=args.orbit_margin)
    o1 = t1 + timedelta(seconds=args.orbit_margin)

    with open(os.path.join(args.out, "orbits.dat"), "w") as fd:
        fd.write("# orbital samples: GPS/Galileo GRG SP3 (15'), BeiDou broadcast (15')\n")
        fd.write(f"# t0 {t0:%Y-%m-%dT%H:%M:%S} GPST\n")
        fd.write("# columns: dt_s sv x_km y_km z_km\n")
        n = 0
        for sv in sorted(used_sv):
            if sv[0] == "C":
                t = o0
                while t <= o1:
                    pos, _, _ = bds_state(bds_select(ephs[sv], t), sv, t)
                    fd.write(f"{(t - t0).total_seconds():g} {sv} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
                    n += 1
                    t += timedelta(seconds=900)
            else:
                for t, pos in sorted(sp3.get(sv, {}).items()):
                    if o0 <= t <= o1:
                        fd.write(f"{(t - t0).total_seconds():g} {sv} {pos[0]:.6f} {pos[1]:.6f} {pos[2]:.6f}\n")
                        n += 1
        print(f"orbits.dat: {n} samples")

    with open(os.path.join(args.out, "clocks.dat"), "w") as fd:
        fd.write("# satellite clocks at observation epochs: GPS/Galileo GRG CLK (30s), BeiDou broadcast\n")
        fd.write("# relativistic correction (-2 r.v/c², or F.e.sqrt(A).sin(E)) is NOT included in clock_s\n")
        fd.write("# group delay: value to apply to the primary signal (BeiDou B1I: -TGD1)\n")
        fd.write(f"# t0 {t0:%Y-%m-%dT%H:%M:%S} GPST\n")
        fd.write("# columns: dt_s sv clock_s relativistic_s group_delay_s\n")
        n = 0
        for (sv, t), (c_s, rel, gd) in sorted(products.items(), key=lambda kv: (kv[0][1], kv[0][0])):
            if t not in epochs:
                continue
            fd.write(f"{(t - t0).total_seconds():g} {sv} {c_s:.12e} {rel:.6e} {gd:.3e}\n")
            n += 1
        print(f"clocks.dat: {n} records")


if __name__ == "__main__":
    main()
