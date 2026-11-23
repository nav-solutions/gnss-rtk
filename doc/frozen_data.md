Frozen real data
================

`src/tests/frozen` replays real observations through the public `Solver` API
and checks performance budgets (availability, horizontal / vertical RMS, max 3D
error) for every technique and constellation.

Data
----

- Rover ESBC00DNK, base MOJN00DNK (65 km baseline), 2020-06-25, from
  [nav-solutions/data](https://github.com/nav-solutions/data):
  `CRNX/V3/{ESBC,MOJN}00DNK_R_20201770000_01D_30S_MO.crx.gz`.
- Signals (match the GRG clock datum): GPS C1W-L1C / C2W-L2W, Galileo E1 / E5a,
  BeiDou B1I / B3I.
- Orbits: GRG SP3 (days 176 & 177, 15') for GPS and Galileo. The GRG SP3 does
  not contain BeiDou: BeiDou orbits are computed from the broadcast ephemerides
  (`NAV/V3/ESBC00DNK_R_20201770000_01D_MN.rnx.gz`, BDS ICD, GEO included) and
  sampled every 15'.
- Clocks: GRG CLK (30 s) for GPS and Galileo, broadcast for BeiDou, given at the
  observation epochs, with the relativistic correction (-2 r.v/c², or F e √A sin E)
  and the group delay to apply (BeiDou B1I: -TGD1).
- Observations without precise products are discarded.
- Truth: antenna reference point (RINEX marker + antenna height).

Two datasets:

| dataset | span | sampling | rover records | size |
|---------|------|----------|---------------|------|
| `dense` | 00:00 - 01:00 | 30 s | 3268 | 0.8 MB |
| `day`   | 00:00 - 23:30 | 5'   | 7221 | 2.0 MB |

At run time, orbits are interpolated (Lagrange, order 10) at the instant the
solver requests (transmission time). Clocks are served at the sampling epoch.

Regenerating / extending
------------------------

```bash
# decompress the observations with crx2rnx (nav-solutions/crx2rnx)
crx2rnx data/CRNX/V3/ESBC00DNK_R_20201770000_01D_30S_MO.crx.gz -o ESBC.rnx
crx2rnx data/CRNX/V3/MOJN00DNK_R_20201770000_01D_30S_MO.crx.gz -o MOJN.rnx

python3 scripts/extract_frozen_data.py \
    --obs ESBC.rnx --obs MOJN.rnx \
    --sp3 data/SP3/C/GRG0MGXFIN_20201760000_01D_15M_ORB.SP3.gz \
    --sp3 data/SP3/C/GRG0MGXFIN_20201770000_01D_15M_ORB.SP3.gz \
    --clk data/CLK/V3/GRG0MGXFIN_20201770000_01D_30S_CLK.CLK.gz \
    --nav data/NAV/V3/ESBC00DNK_R_20201770000_01D_MN.rnx.gz \
    --start "2020-06-25 00:00:00" --duration 3600 --interval 30 \
    --out src/tests/frozen/data/dense

# same with --duration 84600 --interval 300 --out src/tests/frozen/data/day
```

The script only needs the Python standard library. `--systems` restricts the
constellations (default `GEC`).

Budgets and report
------------------

`src/tests/frozen/suite/budgets.rs` holds one budget per case, set to
max(1.4 x, +5 cm) the observed values: any degradation fails, tighten them when
the solver improves.

```bash
cargo test --release frozen
# markdown table (CI job summary)
cargo test --release frozen_performance_report -- --ignored --nocapture
```

Targets (`#[ignore]`): BeiDou CPP (TGD applied to the iono free combination),
multi constellation RTK (single pivot), SPP with an ionosphere model.
