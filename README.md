# SatNOGS Link-Budget Validation Pipeline (MATLAB)

A MATLAB pipeline that **predicts** the received SNR over a satellite pass from
orbital geometry and a station link budget, **measures** the SNR actually seen
in a SatNOGS waterfall recording, and **compares** the two with an offset-fit
method. Target environment: MATLAB R2022b with the Aerospace Toolbox (or
Satellite Communications Toolbox); HDF5 files are read with base MATLAB.

---

## 1. The idea in one page

A ground station receives a satellite signal whose strength is governed by a
handful of physical terms — the **link budget**:

```
Prx = EIRP − FSPL − L_atm − L_line − L_pol + G(el)        [dBW]
N   = 10·log10(k · T_sys · B)                             [dBW]
SNR = Prx − N                                             [dB]
```

| Term | Meaning | Depends on |
|---|---|---|
| `EIRP` | satellite transmit power + antenna gain | the satellite |
| `FSPL` | free-space path loss `20·log10(4π·d·f/c)` | slant range `d` → **time** |
| `L_atm` | atmospheric absorption | elevation → **time** |
| `L_line`, `L_pol` | feedline and polarisation losses | the station |
| `G(el)` | receive antenna gain toward the satellite | elevation → **time** |
| `N` | thermal noise floor, `k·T_sys·B` | the station |

During a pass the slant range swings from ~2500 km at the horizon to a few
hundred km at closest approach, so FSPL — and therefore SNR — traces a
characteristic hill-shaped curve that peaks near maximum elevation. That curve
is *predictable* from a TLE and the station parameters. If the predicted curve
matches what the station actually recorded, the link budget (and your G/T
estimate) is validated.

**Why SNR and never absolute power?** SatNOGS waterfalls are not calibrated.
The HDF5 artifacts store power *rescaled in standard deviations*; the raw
client `.dat` files store dB relative to an arbitrary SDR reference. An
unknown constant is added to every pixel — but it is the *same* constant for
signal and noise bins, so it cancels in the difference `S − N`. SNR is the
only quantity these recordings can support, which is why the whole pipeline
works in SNR.

**Why an offset fit?** Even in SNR, one scalar typically remains unknown or
uncertain (satellite EIRP is rarely published precisely; G/T estimates are
rough). The comparison therefore fits a single constant offset:

- **mean offset** = mean(measured − predicted) — the "fitted station
  constant". It absorbs every unknown scalar in the budget.
- **RMS residual** after removing that offset — how well the *shape* of the
  curve matches. This is the actual validation result.
- **Pearson correlation** — a scale-free second opinion on shape agreement.

Interpretation: offset ≈ 0 with a small, structureless residual means the
budget is validated *absolutely* (your EIRP and G/T numbers are right).
A nonzero offset with a small residual means the *physics* (geometry, FSPL,
elevation dependence) is validated and the offset is your measured correction
to the assumed EIRP/G-over-T. A large or elevation-correlated residual points
at a modelling error (antenna pattern, obstruction, wrong TLE, misalignment).

---

## 2. Repository layout

```
/src                        all code
  run_validation.m          ← the only file you edit to run things
  dbHelpers.m               db/undb conversion handles (single definition)
  generateSyntheticWaterfall.m   synthetic test waterfall (selftest)
  extractMeasuredSNR.m      measured SNR from a waterfall
  compareAndPlot.m          offset fit, statistics, figures
  computePassGeometry.m     TLE → az/el/range/Doppler histories
  predictSNR.m              link budget → predicted SNR(t)
  loadStationConfig.m       reads/validates the station config
  parseWaterfall.m          format dispatcher (.h5/.hdf5/.dat)
  parseWaterfallHDF5.m      SatNOGS .h5 artifacts (dequantises uint8 data)
  parseWaterfallDAT.m       (step 4, not yet built) raw client .dat files
/config
  station_uhf.m             your station parameters — EDIT THIS
/data                       drop your .h5 / .dat / .tle files here
/output                     figures (.png + .fig) and results (.mat) land here
```

**Conventions used everywhere:** power and SNR in dB; frequency in Hz; time as
`datetime` in UTC; angles in degrees. The helpers `db(x) = 10·log10(x)` and
`undb(x) = 10^(x/10)` are defined once in `src/dbHelpers.m`.

---

## 3. Quick start

Everything is driven from `src/run_validation.m`. Open it, edit the config
block at the top, press **F5**. The `MODE` variable selects what runs:

| MODE | What it does | Needs |
|---|---|---|
| `"selftest"` | end-to-end test of the measurement chain on synthetic data | nothing — no toolboxes, no files |
| `"geomtest"` | pass table + predicted SNR from a real TLE | Aerospace Toolbox, a TLE in `/data`, edited config |
| `"hdf5"` | full pipeline on a SatNOGS network artifact | *(build step 3/5 — not yet available)* |
| `"dat"` | full pipeline on a raw station `.dat` file | *(build step 4 — not yet available)* |

### 3.1 First run: the selftest

```matlab
MODE = "selftest";   % in run_validation.m, then F5
```

No data or toolboxes needed. It fabricates a 10-minute overhead pass
(elevation −5° → 60° → −5°, Doppler +10 kHz → −10 kHz), builds a synthetic
waterfall containing a Doppler-shifted signal ridge with a *known* SNR
profile, then runs the real extraction and comparison code against it.

It **passes** when the extracted SNR matches the known truth within 0.5 dB
(max error above 5° elevation) and correlation exceeds 0.99. Expect roughly:

```
Selftest: max |meas - true| = 0.196 dB (limit 0.50), correlation = 0.9999
SELFTEST PASSED.
```

Run this first, and re-run it after any code change — it proves the SNR logic
with zero external dependencies. Three diagnostic figures land in `/output`;
figure 3 (waterfall with the red Doppler track on the ridge) is the visual
proof that time/frequency alignment works.

### 3.2 Configure your station

Edit `config/station_uhf.m`. Every `EDIT ME` field matters:

```matlab
station.lat_deg          = 52.0;   % station latitude  (deg, +N)
station.lon_deg          = 4.4;    % station longitude (deg, +E)
station.alt_m            = 10;     % altitude (m)
station.freq_Hz          = 435e6;  % downlink frequency (Hz)
station.ant_gain_dBi     = 12;     % antenna boresight gain (dBi)
station.ant_pattern      = [];     % optional @(el_deg) -> dBi
station.sys_noise_temp_K = 500;    % system noise temperature (K)
station.rx_bw_Hz         = 5000;   % SNR noise bandwidth (Hz)
station.line_loss_dB     = 1;      % feedline + connectors (dB)
station.pol_loss_dB      = 3;      % polarisation mismatch (dB)
station.sat_eirp_dBW     = -3;     % satellite EIRP (dBW); 0.5 W ≈ -3 dBW
```

Guidance:

- **`sys_noise_temp_K`** — dominated by the LNA and sky noise at UHF; 400–800 K
  is typical for a SatNOGS-class station without careful optimisation.
- **`rx_bw_Hz`** — keep equal to the `sig_bw_Hz` used in the extraction
  (default 5000 Hz) so predicted and measured SNR share the same bandwidth.
  This matters: SNR is bandwidth-dependent.
- **`ant_pattern`** — leave `[]` for constant gain, or supply e.g.
  `@(el) 12 - 0.002*(90-el).^2` for a crude elevation roll-off.
- **`sat_eirp_dBW`** — your best guess. If it's wrong by X dB, the comparison
  simply reports a fitted offset of −X dB; the shape validation is unaffected.

**Precedence rule:** when the waterfall is a SatNOGS HDF5 artifact, its
embedded metadata (station location, frequency, TLE) *overrides* the config
for the geometry — the artifact knows exactly where and what was recorded.
The config always supplies the RF terms (gains, losses, T_sys, EIRP), which
the artifact does not carry.

### 3.3 Check the geometry: geomtest

1. Get a **current** TLE for your satellite (Celestrak, or the SatNOGS
   observation page — TLEs go stale within days). Save the 3 lines (name +
   line 1 + line 2) as `data/target.tle`.
2. In `run_validation.m` set `MODE = "geomtest"` and set `obs_start` /
   `obs_stop` to a window that contains a pass over your station.
3. F5. You get a pass table every 30 s:

```
          time (UTC)   az deg   el deg   range km  doppler kHz     SNR dB
 2026-01-01 12:03:30    201.3      4.2     1890.4         9.61       18.9
 2026-01-01 12:05:30    247.8     34.7      750.2         3.12       26.8
 ...
```

Sanity checks: elevations should match any pass predictor (gpredict, n2yo,
the SatNOGS page itself) to a fraction of a degree; Doppler should swing
roughly ±10 kHz at 435 MHz, crossing zero at maximum elevation; SNR should
peak at max elevation. If max elevation is negative you get a clear error —
wrong TLE or wrong time window.

### 3.4 Real observations (hdf5 / dat modes)

The HDF5 parser is built (step 3): `wf = parseWaterfall('data/....h5')`
returns the unified waterfall struct, with the artifact's TLE / frequency /
station location in `wf.meta`. Note that artifacts store the spectrum as
uint8 with per-channel offset/scale vectors; the parser dequantises back to
uncalibrated dB automatically (per-channel offsets would NOT cancel in the
S − N difference, so this step is mandatory — see the header comment).

The end-to-end modes arrive with steps 4–5: download the waterfall artifact
(.h5) from an observation page on network.satnogs.org into `/data`, point
`paths.wf_file` at it, set `MODE = "hdf5"`, F5. For raw `.dat` files from
your own station (`/tmp/.satnogs/data/` on the Pi) you additionally supply
the TLE and pass window yourself, since `.dat` files carry no metadata.

---

## 4. How each module works (education section)

### Prediction side

**`computePassGeometry`** builds a `satelliteScenario`, propagates the TLE
(SGP4 under the hood), and asks for azimuth/elevation/range histories with
`aer(groundStation, satellite)`. Doppler is *not* read from the toolbox — it
is computed from the numerical range rate:

```
rr = gradient(range) / gradient(t)          % m/s, + = receding
doppler = −rr/c · f                         % Hz, + = approaching (rising pass)
```

This is the classic non-relativistic Doppler formula; at LEO velocities
(|rr| ≤ ~7 km/s) that is exact to well under 1 Hz.

**`predictSNR`** evaluates the link-budget equation from §1 per time step.
Details worth knowing:

- FSPL uses the exact slant range from the geometry, so the elevation
  dependence of the prediction comes almost entirely from FSPL.
- Atmospheric loss uses `gaspl` (ITU-R P.676 gas model) with a standard
  atmosphere and a cosecant path **only above 1 GHz** — that is gaspl's
  validity floor. At 435 MHz atmospheric absorption is genuinely tiny, so
  the code falls back to a fixed 0.2 dB and tells you so with a warning.
- Samples below 5° elevation are flagged `excluded`: near the horizon,
  multipath, obstructions and refraction make both prediction and
  measurement unreliable.
- The constant terms are returned in a `budget` struct — that is your thesis
  link-budget table.

### Measurement side

**`extractMeasuredSNR`** walks the waterfall row by row (each row = one
spectrum at one time):

1. **Align time.** Waterfall rows are seconds-since-start; the geometry is
   absolute UTC. If the waterfall has a `start_time` (HDF5 artifacts do), the
   two are anchored exactly; otherwise both are assumed to start together.
   A `time_offset_s` parameter lets you nudge the alignment manually — watch
   figure 3 to see when it's right.
2. **Follow the Doppler track.** The predicted Doppler is interpolated onto
   each row time; the signal window (default ±2.5 kHz) is centred on it. This
   is why measurement needs geometry: the signal moves through the passband
   during the pass, and a fixed window would slide off it.
3. **Signal estimate `S`** = mean of the top-3 bins in the window (robust to
   the exact carrier position within the window).
4. **Noise estimate `N`** = *median* of the guard band 20–50 kHz either side
   of the track — far enough to be clean of signal, and the median is immune
   to interference spikes in a way the mean is not.
5. `SNR = S − N`. The unknown calibration constant cancels here (§1).

**`compareAndPlot`** interpolates the prediction onto the measurement
timestamps, masks to elevation > 5°, and computes the offset / RMS residual /
correlation described in §1. It saves three figures (PNG + .fig) and a
results `.mat` to `/output`:

1. **SNR vs time** — predicted (solid) and measured (dots), elevation on the
   right axis. The headline plot.
2. **Residuals vs elevation** — after offset removal. Should be a
   structureless cloud around zero; a trend with elevation means the antenna
   pattern or horizon model is wrong.
3. **Waterfall + Doppler track** — the red line must lie on the signal ridge.
   If it is parallel but shifted in time, adjust `time_offset_s`; if shifted
   in frequency, the TLE is stale or the transmitter is off-frequency.

### Test harness

**`generateSyntheticWaterfall`** builds a fake waterfall with a known answer:
Gaussian noise floor (in dB — real rows are averages of many FFTs, which
makes the noise statistics near-Gaussian in dB) plus a 3-bin signal ridge at
`center + doppler(t)` whose height is set so that the extraction definition
`S − N` recovers the target SNR exactly. Feeding a known profile through the
*real* extraction and comparison code and demanding ≤ 0.5 dB error is the
acceptance test that protects every later change.

---

## 5. Reading the results

Console output of a real run ends with a paragraph like:

> 412 samples above 5 deg elevation. Mean offset (meas−pred) = −6.31 dB;
> RMS residual after offset removal = 1.42 dB; Pearson correlation = 0.964.

How to report this:

- **Offset −6.3 dB** → your assumed EIRP·G/T product is 6.3 dB optimistic.
  Quote it as the fitted station constant, or fold it into an implied G/T:
  `G/T_implied = G/T_assumed + offset`.
- **RMS 1.4 dB** → the shape of the SNR curve is reproduced to ~1.4 dB — for
  an uncalibrated crowdsourced receiver that is a strong validation.
- **Correlation 0.96** → prediction and measurement rise and fall together.

Everything needed to regenerate the numbers is saved in
`output/<label>_results.mat` (per-sample vectors, mask, residuals).

---

## 6. Troubleshooting

| Symptom | Likely cause / fix |
|---|---|
| `satelliteScenario not found` | Aerospace / SatCom Toolbox missing — geomtest and real modes need it; selftest does not. |
| `Satellite never rises above the horizon` | Wrong TLE, wrong station coordinates, or the obs window doesn't contain the pass. |
| Doppler track (fig. 3) misses the ridge in **time** | Adjust `time_offset_s` in `run_validation.m`. |
| Track misses in **frequency** | Stale TLE (get one from the observation date) or transmitter offset. |
| `gaspl unavailable or f < 1 GHz` warning | Expected at UHF — 0.2 dB fixed loss is used; harmless. |
| Noise-guard fallback warning | Waterfall span too narrow for the 20–50 kHz guard; the code falls back to all bins away from the track. Fine, but consider narrower guards. |
| Selftest fails after an edit | You changed the extraction/comparison logic — fix before trusting real results. |

---

## 7. Build status

| Step | Content | Status |
|---|---|---|
| 1 | selftest chain (synthetic waterfall → extraction → comparison) | ✅ merged |
| 2 | pass geometry + link-budget prediction (geomtest) | ✅ merged |
| 3 | SatNOGS HDF5 artifact parser + dispatcher | ✅ this PR |
| 4 | raw client `.dat` parser | ⏳ |
| 5 | full integration in hdf5 mode | ⏳ |
