# FX4 fly scan — uncertainties and the negative-intensity lift

Status: **implemented in Matilda (2026-11)**; Igor not yet. Porting steps for
Igor: `fx4-reduction-guide.md` §3.7, §3.8, §3.11, §6.3. Code:
`matilda/fx4support.py` (constants, `current_sigma`), `supportFunctions.py`
(`_calculatePDErrorFlyFX4`, `smooth_r_data`, `rebin_log_data`,
`calibrateAndSubtractFlyscan`). Tests: `tests/test_fx4_noise.py`.

## 1. Test data

`TestData/error estimates/` (not committed) and `TestData/FX4Set/flyscan/`:

| pair | energy | character |
|---|---|---|
| `GC_SRM3600_0036` / `AirBlank_0035` | 12 keV | weak: 2–7 pA UPD over most of Q |
| `C2p5Top_0010` / `NMRTubeBlank_0003` | 24 keV | very strong (4e5 pA at Q = 5e-4), MSAXS = 7.0, T = 0.38 |
| `GC_SRM3600_0131` / `Blank_0130` | 21 keV | commissioning pair, 20–120 pA |
| `AirBlank_28keV_220_0008` | 28 keV | blank only; recorded dark current 4.53 pA is **above** the 2.22 pA high-Q signal (see §5) |

## 2. How it is measured

**Halves test.** Reduce the even and the odd raw points independently, compare
ln I of the two SMR curves (half B on half A's Q by a log–log cubic spline),
and divide the robust spread of the difference by the combined quoted error:
**true/quoted**, target 1. Both halves use the *full* scan's beam centre and
peak maximum, otherwise each half fits its own peak and that ~1 % shows up as
fake noise. This is `tests/test_fx4_noise.py::_halves_ratio`.

Do not judge error bars by the visible scatter of the SMR curve: smoothed
neighbours share raw data, so correct error bars look 2–5× too large to the
eye. Raw-level checks use second differences of the unsmoothed signal, where
neighbouring points are uncorrelated (lag-1 autocorrelation ≈ 0).

## 3. Results

True/quoted at the SMR stage, per Q band:

| Q range | GC_0036 before | C2p5 before | GC_0131 before | GC_0036 now | C2p5 now | GC_0131 now |
|---|---|---|---|---|---|---|
| 3e-4 – 1e-3 | – | 21 | – | – | 1.00 | – |
| 1e-3 – 3e-3 | 0.27 | 11 | 1.17 | 1.07 | 0.84 | 1.05 |
| 3e-3 – 1e-2 | 0.14 | 6.3 | 0.93 | 0.75 | 0.93 | 0.82 |
| 1e-2 – 3e-2 | 0.22 | 4.2 | 0.59 | 1.12 | 0.86 | 0.72 |
| 3e-2 – 0.1 | 0.18 | 4.8 | 0.93 | 0.75 | 1.60 | 0.96 |
| 0.1 – 0.3 | 0.18 | 4.7 | 0.69 | 0.76 | 1.32 | 0.77 |

Before, errors were ~5× too large on the weak sample and 5–20× too small on
the strong one. The commissioning pair looked right only because two errors
cancelled at its 20–100 pA signal level.

## 4. What was wrong, and what the code does now

### 4.1 Per-interval error (FX4)

Was `upd_sigma/sqrt(N)`, N = reads per interval. The raw noise follows

```
sigma^2 = upd_sigma^2/(N·values_per_read)  +  q_eff·max(I − dark, 0)/t
```

* **Floor.** Overstated by ~√10 before.
* **Shot-like term**, missing before, dominant from ~1 pA net. Fitted
  `q_eff`: 6.0–7.5e-4 pC at 12 keV, 6.4e-4 at 21, 6.5–9e-4 at 24, 4.7–8.9e-4 at
  28 keV. E/3.65 eV would give 5.3e-4 → 1.2e-3 over that range, so the fitted
  value does not follow it; one constant, **6.5e-4 pC**, is used. Diode
  efficiency does not enter (the term is in measured current).
* On the strong sample at 3e3–3e4 pA the raw relative noise is 0.17–0.20 %,
  fully explained by the two terms. About **0.1 % of I** remains after all
  averaging at the lowest Q, added once at the end.

All terms depend on signal level only, never on Q, so they carry over to
samples 10× weaker or 1000× stronger than glassy carbon.

### 4.2 Smoothing and rebinning

`smooth_r_data` averages w raw points (41–47 below 30 pA, 7–15 above, at
t ≈ 0.019 s) and left the error unchanged; now it divides by √w and records w.
Rebinning then counts only independent points, RMS/√max(n/w, 1). The old
"error of the mean" switch divided by √n and double-counted.

Smoothing is now done in linear intensity for FX4. Averaging ln(I) read 11–12 %
low at high Q on C2p5 and 1–2 % low on glassy carbon; on strong, steep data
linear and log differ by < 0.07 %.

### 4.3 `SMR_Error *= PeakToPeakTransmission` (both chains)

Shrank the relative error by T_pp: 0.92–0.94 on ordinary samples, 0.055 on
C2p5 (18×). It came from Igor, where it cancels an extra `/T` in Igor's own
formula; Matilda has no such extra division. Removed.

## 5. The negative-intensity lift

The lift (Igor `IN3_FixNegativeIntensities`) adds one constant to the whole
curve when the high-Q intensity goes negative, so that log(I) works. A constant
is absorbed by the flat background of every model used downstream, so it does
not bias fits. What was wrong was **where** it ran on FX4:

* It ran on raw data. Weak FX4 points scatter around zero, so the most negative
  point was a noise spike: lifts of 0.3–0.7 pA against a net high-Q signal of
  0.3–0.5 pA, different for sample and blank and from scan to scan.
* The log-domain smoothing that followed then averaged ln(I) of data with a
  random offset and noise comparable to the signal.

Now the FX4 chain smooths first (linear) and lifts after, only if the smoothed
curve still goes negative. Then the lift follows the background level, not a
spike. On `GC_SRM3600_0131` this is a flat +0.17 cm⁻¹ at every Q.

The 28 keV air blank is the case the lift exists for: its recorded dark
current (4.53 pA) is above the measured high-Q signal (2.22 pA), so every
high-Q point is about −2.3 pA after subtraction. That is a dark-current
measurement problem at the beamline, worth checking; the lift handles it.

## 6. Open items

* **Below 12 keV:** no data yet. q_eff holds 12–28 keV.
* **What `upd_sigma` is:** confirm against `bits_usaxs/docs/FX4_data_formats.md`.
  The `N·values_per_read` floor is inferred from data.
* **Igor:** port per `fx4-reduction-guide.md` §3.11; also revisit Igor's own
  `/3` (12/2013) and the SMR_Error not divided by MSAXSCorrection.
* **Scaler chain:** it still averages ln(I) and lifts before smoothing. The same
  weak-signal bias probably applies there; not changed, to keep historical
  results.
* **Dark current at 28 keV:** `AirBlank_28keV_220_0008` records a dark current
  2.3 pA above its own high-Q signal; worth checking how it was measured.
