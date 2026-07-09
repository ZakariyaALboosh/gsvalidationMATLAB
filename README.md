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
  undoDopplerCorrection.m   utility: corrected waterfall -> raw sky view
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
| `"hdf5"` | full pipeline on a SatNOGS network artifact | Aerospace Toolbox, a .h5 in `/data`, edited config |
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

**hdf5 mode is the full pipeline**: download the waterfall artifact (.h5)
from an observation page on network.satnogs.org into `/data`, point
`paths.wf_file` at it, set `MODE = "hdf5"`, F5. Everything the geometry
needs (TLE, frequency, station location, start time) comes from the
artifact's own metadata and overrides the config; the config supplies the
RF terms. The parser dequantises the uint8 spectrum back to uncalibrated dB
automatically (artifacts store per-channel offset/scale vectors, which
would NOT cancel in the S − N difference — see the parser header).

**Crucial: SatNOGS waterfalls are Doppler-corrected.** The client shifts
the receive frequency along the predicted Doppler curve while recording, so
in the waterfall the *satellite* appears as a near-vertical trace at a
constant offset (transmitter off-frequency + oscillator drift), while
*fixed terrestrial carriers* (birdies) appear as inverted S-curves — the
mirror image of the correction. This was verified empirically on a real
artifact: a strong carrier followed −doppler(t) to within one bin. Hence
`wf_doppler_corrected = true` in the config block (the extraction tracks a
constant offset); `freq_offset_Hz` nudges that offset if the transmitter is
off-frequency — tune it by eye on the waterfall/track overlay figure.

**If your file is NOT Doppler-corrected** (a non-SatNOGS recording, or an
unusual client build): set `wf_doppler_corrected = false` in the config
block — the extraction then follows the full ±10 kHz Doppler S-curve
instead of a constant offset. Both cases are first-class; the setting just
picks the track model. To *check* which kind you have, look at the
waterfall figure: satellite vertical + birdies S-curved ⇒ corrected;
satellite S-curved ⇒ uncorrected.

**Converting a corrected waterfall back to the raw sky view**: the utility
`undoDopplerCorrection` shifts every row by the predicted Doppler, so the
satellite becomes the classic S-curve again and birdies become vertical
(verified on the real artifact: the birdie's residual drops to 30 Hz,
under one bin). Optionally writes the result as an .h5 readable by this
pipeline:

```matlab
wf  = parseWaterfall('data/observation.h5');
geo = computePassGeometry(wf.meta.tle, wf.meta.lat, wf.meta.lon, ...
          wf.meta.alt, wf.start_time, ...
          wf.start_time + seconds(wf.t_s(end)), 1, wf.meta.frequency_Hz);
wf_raw = undoDopplerCorrection(wf, geo, ...
          struct('out_h5', 'output/observation_raw.h5'));
```

Bins that shift in from outside the recorded passband are filled with the
row median — treat the outer ±|doppler| edge strips as synthetic.

Pick a **well-vetted, strong observation** for validation runs: on a weak
or QRM-dominated recording the measured "SNR" sits at the estimator's noise
bias floor (top-3-of-window on pure noise gives a few dB) and correlation
with the prediction collapses — the figures make this immediately obvious.

The `.dat` mode arrives with step 4: for raw files from your own station
(`/tmp/.satnogs/data/` on the Pi) you additionally supply the TLE and pass
window yourself, since `.dat` files carry no metadata.

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
2. **Follow the signal track.** The track model is
   `track(t) = doppler_factor · doppler(t) + freq_offset_Hz`. For a raw
   uncorrected spectrum `doppler_factor = 1`: the signal sweeps ±10 kHz
   through the passband and the window must follow it. For SatNOGS
   waterfalls (Doppler-corrected, see §3.4) `doppler_factor = 0`: the
   satellite sits at a near-constant offset. The signal window (default
   ±2.5 kHz) is centred on the track either way.
3. **Signal estimate `S`** = mean of the top-3 bins in the window (robust to
   the exact carrier position within the window).
4. **Noise estimate `N`** = *median* of the guard band 20–50 kHz either side
   of the track, **restricted to the inner 80 % of the passband** — real SDR
   waterfalls roll off 10–30 dB at the band edges (anti-alias filter), and
   guard bins out there would fake a low noise floor (+16 dB SNR bias seen
   on a real artifact). If fewer than 64 usable guard bins remain, the
   estimator falls back to the median of the central 60 % of the passband,
   clear of the signal window. The median is immune to the narrow signal
   and to interference spikes in a way the mean is not.
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
| Track misses in **frequency** | For hdf5: transmitter off-frequency — set `freq_offset_Hz`. For raw spectra: stale TLE (get one from the observation date). |
| Strong **inverted S-curve** in the waterfall | A fixed terrestrial carrier mirrored by the Doppler correction — not the satellite. The satellite is the near-vertical trace. |
| Measured SNR flat ≈ few dB, correlation ≈ 0 | No usable satellite signal on the track (weak/QRM-dominated observation): you are seeing the top-3-of-noise bias floor. Pick a stronger, well-vetted observation. |
| `gaspl unavailable or f < 1 GHz` warning | Expected at UHF — 0.2 dB fixed loss is used; harmless. |
| Noise-guard fallback warning | The 20–50 kHz guard falls outside the flat part of the passband; the code uses the central-passband median instead. Expected on 48 kHz-wide waterfalls; harmless. |
| Selftest fails after an edit | You changed the extraction/comparison logic — fix before trusting real results. |

---

## 7. Build status

| Step | Content | Status |
|---|---|---|
| 1 | selftest chain (synthetic waterfall → extraction → comparison) | ✅ merged |
| 2 | pass geometry + link-budget prediction (geomtest) | ✅ merged |
| 3 | SatNOGS HDF5 artifact parser + dispatcher | ✅ merged |
| 4 | raw client `.dat` parser | ⏳ next |
| 5 | full integration in hdf5 mode | ✅ this PR |
