"""Independent equation checks, NOT execution of MATLAB source or patient data."""
from pathlib import Path
import json
import re
import numpy as np
from scipy.signal import fftconvolve
from scipy.signal.windows import dpss

ROOT = Path(__file__).resolve().parents[1]
checks = {}


def fir(fs, band, seconds):
    order = max(4, 2 * int(np.ceil(seconds * fs / 2)))
    n = np.arange(-order // 2, order // 2 + 1)
    lo, hi = band
    h = 2 * hi / fs * np.sinc(2 * hi / fs * n)
    if lo > 0:
        h -= 2 * lo / fs * np.sinc(2 * lo / fs * n)
    h *= 0.54 - 0.46 * np.cos(2 * np.pi * np.arange(order + 1) / order)
    gain = h.sum() if lo == 0 else np.sum(h * np.cos(2 * np.pi * np.sqrt(lo * hi) / fs * n))
    return h / gain


def apply_fir(x, h):
    valid = np.isfinite(x)
    y = fftconvolve(np.where(valid, x, 0), h, mode="same")
    support = fftconvolve(valid.astype(float), np.ones(len(h)), mode="same")
    ok = support >= len(h) - 0.01
    y[~ok] = np.nan
    return y, ok


def spectrum(y, fs, seconds=4):
    n = round(seconds * fs)
    nfft = 2 ** int(np.ceil(np.log2(n)))
    w = 0.5 - 0.5 * np.cos(2 * np.pi * np.arange(n) / n)
    w /= np.linalg.norm(w)
    reps = []
    for a in range(0, len(y) - n + 1, n // 2):
        v = y[a:a + n]
        if not np.isfinite(v).all():
            continue
        z = np.fft.rfft((v - v.mean()) * w, nfft) / np.sqrt(fs)
        z[1:-1] *= np.sqrt(2)
        reps.append(z)
    return np.mean(np.abs(reps) ** 2, axis=0), np.fft.rfftfreq(nfft, 1 / fs)


fs = 100
t = np.arange(400) / fs
p, f = spectrum(np.cos(2 * np.pi * 6 * t), fs)
area = np.trapezoid(p, f)
assert abs(area - 0.5) < 1e-6
checks["one_sided_PSD_area"] = float(area)

n, nw, k = 100, 1.5, 2
ix = np.arange(n)
diag = ((n - 1 - 2 * ix) / 2) ** 2 * np.cos(2 * np.pi * nw / n)
off = np.arange(1, n) * (n - np.arange(1, n)) / 2
T = np.diag(diag) + np.diag(off, 1) + np.diag(off, -1)
eig, V = np.linalg.eigh(T)
W = V[:, np.argsort(eig)[::-1][:k]]
reference = dpss(n, nw, Kmax=k).T
error = np.max(np.abs(np.abs(W.T @ reference) - np.eye(k)))
assert error < 1e-10
checks["DPSS_match_to_SciPy_max_error"] = float(error)

fs = 1000
t = np.arange(12000) / fs
x = (1 + 0.65 * np.cos(2 * np.pi * 6 * t)) * np.cos(2 * np.pi * 90 * t)
hB = fir(fs, (20, 450), .15)
hL = fir(fs, (0, 40), .5)
xb, _ = apply_fir(x, hB)
e, ok = apply_fir(np.abs(xb), hL)
yb, _ = apply_fir(2.7 * x, hB)
ey, oky = apply_fir(np.abs(yb), hL)
scaling_error = np.max(np.abs(ey[ok & oky] - 2.7 * e[ok & oky]))
assert scaling_error < 1e-10
checks["raw_to_envelope_scaling_max_error"] = float(scaling_error)
gap = x.copy()
gap[5000] = np.nan
_, valid_gap = apply_fir(gap, hB)
assert not valid_gap[5000 - len(hB) // 2:5001 + len(hB) // 2].any()
checks["gap_support_is_not_concatenated"] = True

ea = np.interp(np.arange(1200) / 100, t, e)
p, f = spectrum(ea, 100)
noise = np.full_like(p, np.nan)
for i, hz in enumerate(f):
    if hz < 2 or hz > 25:
        continue
    use = (np.abs(f - hz) > max(.75, f[1] - f[0])) & (np.abs(f - hz) <= 5) & (f >= 2) & (f <= 25)
    values = p[use & np.isfinite(p) & (p > 0)]
    if len(values) >= 4:
        noise[i] = np.quantile(values, .5)
ratio = p / np.maximum(noise, max(np.nanmax(p) * 1e-12, np.finfo(float).tiny))
pick = np.flatnonzero((f >= 4) & (f <= 12) & np.isfinite(ratio))
frequency = f[pick[np.argmax(ratio[pick])]]
assert abs(frequency - 6) < .5
checks["Welch_CFAR_tone_frequency_Hz"] = float(frequency)

# Masking a target block cannot alter the training observer beyond its support.
changed = x.copy()
changed[(t >= 4) & (t < 7)] += 20 * np.sin(2 * np.pi * 31 * t[(t >= 4) & (t < 7)])
z, _ = apply_fir(changed, hB)
q1, qok1 = apply_fir(xb ** 2, hL)
q2, qok2 = apply_fir(z ** 2, hL)
train = qok1 & qok2 & ((t < 3.4) | (t > 7.6))
leakage = np.max(np.abs(q1[train] - q2[train]))
assert leakage < 1e-10
checks["fixed_observer_guard_leakage_max_error"] = float(leakage)

rng = np.random.default_rng(42)
y = rng.normal(size=1000)
s = .6 * y + rng.normal(size=1000)
r = y - s
identity = np.var(y) - np.var(s) - np.var(r) - 2 * np.mean((s - s.mean()) * (r - r.mean()))
assert abs(identity) < 1e-12
checks["variance_cross_term_identity_error"] = float(abs(identity))

t = np.arange(3000) / 100
g = ((t >= 2) & (t < 4)) | ((t >= 12) & (t < 14)) | ((t >= 22) & (t < 24))
s = 10 * g * np.cos(2 * np.pi * 6 * t)
assert abs(np.mean(s ** 2) - g.mean() * np.mean(s[g] ** 2)) < 1e-12
checks["occupancy_power_identity"] = True

# The mean-variance formula behind the Gaussian AM calibration, by Monte Carlo.
fs = 1000
t = np.arange(6000) / fs
m = 1 + .65 * np.cos(2 * np.pi * 6 * t)
h = fir(fs, (20,450), .15)
expected_var = fftconvolve(m ** 2 + .3 ** 2, h ** 2, mode="same")
sample = m[:, None] * rng.normal(size=(len(t), 600)) + .3 * rng.normal(size=(len(t),600))
filtered = fftconvolve(sample, h[:, None], mode="same", axes=0)
select = np.arange(500,5500)
relative_error = abs(np.mean(filtered[select] ** 2) / np.mean(expected_var[select]) - 1)
assert relative_error < .02
checks["Gaussian_AM_variance_Monte_Carlo_relative_error"] = float(relative_error)

package_files = list((ROOT / "+tbench").glob("*.m"))
known = {p.stem for p in package_files}
missing = set()
for path in ROOT.rglob("*.m"):
    missing |= set(re.findall(r"tbench\.([A-Za-z][A-Za-z0-9_]*)\s*\(", path.read_text())) - known
assert not missing, missing
checks["package_direct_calls_resolve"] = True
checks["matlab_files"] = len(list(ROOT.rglob("*.m")))
checks["MATLAB_native_tests_executed_here"] = False
checks["patient_data_analyzed_here"] = False
(ROOT / "validation" / "numerical_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
print(json.dumps(checks, indent=2))
