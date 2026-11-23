"""Independent numerical reference checks, NOT execution of the MATLAB source."""
from pathlib import Path
import json
import re
import numpy as np
from scipy.signal import fftconvolve
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parents[1]
checks = {"MATLAB_native_tests_executed_here": False, "patient_data_analyzed_here": False}


def design_aa(fs, pass_hz, stop_hz, attenuation=90, ripple_db=.05):
    beta = .1102 * max(0, attenuation - 8.7)
    fc = (pass_hz + stop_hz) / 2
    order = 2 * int(np.ceil(max(8, (attenuation - 8) /
                               (2.285 * 2 * np.pi * (stop_hz - pass_hz) / fs)) / 2))
    for _ in range(12):
        n = np.arange(-order // 2, order // 2 + 1)
        h = 2 * fc / fs * np.sinc(2 * fc / fs * n) * np.kaiser(order + 1, beta)
        h /= h.sum()
        nfft = 2 ** int(np.ceil(np.log2(max(16384, 16 * len(h)))))
        H = np.abs(np.fft.rfft(h, nfft))
        f = np.fft.rfftfreq(nfft, 1 / fs)
        ripple = np.max(np.abs(20 * np.log10(np.maximum(H[f <= pass_hz], np.finfo(float).tiny))))
        rejection = -20 * np.log10(np.max(H[f >= stop_hz]))
        if ripple <= ripple_db and rejection >= attenuation:
            return h, float(ripple), float(rejection)
        order = 2 * int(np.ceil(order * 1.15 / 2))
    raise AssertionError("Filter failed its declared response")


def apply_fir(x, h):
    if x.ndim == 1:
        x = x[:, None]
    finite = np.isfinite(x)
    y = fftconvolve(np.where(finite, x, 0), h[:, None], mode="same", axes=0)
    support = fftconvolve(finite.astype(float), np.ones((len(h), 1)), mode="same", axes=0)
    y[support < len(h) - .01] = np.nan
    return y


def rate_change(x, fs, target, pass_hz, stop_hz):
    h, ripple, rejection = design_aa(fs, pass_hz, stop_hz)
    filtered = apply_fir(np.asarray(x), h)
    t = np.arange(int(np.floor((len(x) - 1) / fs * target)) + 1) / target
    position = t * fs
    if np.all(np.abs(position - np.round(position)) < 1e-7):
        return filtered[np.round(position).astype(int)], t, ripple, rejection
    length = 20
    index = np.floor(position).astype(int)[:, None] + np.arange(-length + 1, length + 1)
    distance = position[:, None] - index
    w = np.sinc(distance) * np.i0(8.6 * np.sqrt(np.maximum(0, 1 - (distance / length) ** 2))) / np.i0(8.6)
    w[np.abs(distance) > length] = 0
    w /= w.sum(axis=1, keepdims=True)
    inside = np.all((index >= 0) & (index < len(x)), axis=1)
    index = np.clip(index, 0, len(x) - 1)
    y = np.full((len(t), filtered.shape[1]), np.nan)
    for m in range(filtered.shape[1]):
        v = filtered[index, m]
        valid = inside & np.all(np.isfinite(v) | (np.abs(w) < 1e-12), axis=1)
        z = np.sum(np.where(np.isfinite(v), v, 0) * w, axis=1)
        y[valid, m] = z[valid]
    return y, t, ripple, rejection


def tone_estimate(y, t, frequency):
    good = np.isfinite(y)
    X = np.column_stack((np.ones(len(t)), np.cos(2 * np.pi * frequency * t),
                         np.sin(2 * np.pi * frequency * t)))
    b = np.linalg.lstsq(X[good], y[good], rcond=None)[0]
    return np.linalg.norm(b[1:]), np.arctan2(-b[2], b[1])


responses, tone_rows = [], []
for fs in (1000, 1111, 2500):
    h, ripple, rejection = design_aa(fs, 30, 45)
    assert np.max(np.abs(h - h[::-1])) < 1e-14
    responses.append({"inputFs": fs, "taps": len(h), "passErrorDB": ripple,
                      "stopAttenuationDB": rejection, "halfSupportSeconds": (len(h) - 1) / (2 * fs)})
    t = np.arange(12 * fs) / fs
    for frequency in (3, 6, 9, 18, 27):
        y, u, _, _ = rate_change(np.cos(2 * np.pi * frequency * t), fs, 100, 30, 45)
        gain, phase = tone_estimate(y[:, 0], u, frequency)
        assert abs(gain - 1) < .006 and abs(phase) < 1e-3
        tone_rows.append({"inputFs": fs, "frequency": frequency, "gain": float(gain), "phase": float(phase)})
    y, u, _, _ = rate_change(np.cos(2 * np.pi * 90 * t), fs, 100, 30, 45)
    gain, _ = tone_estimate(y[:, 0], u, 90)
    assert gain < 5e-5
checks["AA_response"] = responses
checks["passband_tones_and_harmonics"] = tone_rows

fs = 2500
t = np.arange(12 * fs) / fs
x = np.cos(2 * np.pi * 150 * t + (10 / 6) * np.sin(2 * np.pi * 6 * t))
y, u, _, _ = rate_change(x, fs, 1250, 450, 537.5)
expected = np.cos(2 * np.pi * 150 * u + (10 / 6) * np.sin(2 * np.pi * 6 * u))
good = np.isfinite(y[:, 0])
fm_error = np.sqrt(np.mean((y[good, 0] - expected[good]) ** 2))
assert fm_error < 1e-3
checks["pure_FM_carrier_RMS_error"] = float(fm_error)

# Fractional-rate gaps and shared clocks: unsupported intervals stay missing.
fs = 1111
t = np.arange(12 * fs) / fs
x = np.column_stack((np.cos(2 * np.pi * 6 * t), 2 * np.cos(2 * np.pi * 6 * t)))
x[5000] = np.nan
y, u, _, _ = rate_change(x, fs, 100, 30, 45)
assert not np.isfinite(y[np.abs(u - 5000 / fs) < .1]).any()
good = np.isfinite(y).all(axis=1)
assert np.max(np.abs(y[good, 1] - 2 * y[good, 0])) < 1e-12
checks["fractional_rate_gaps_and_channel_clock"] = True

# Demodulate before reducing: a raw carrier cannot be moved to 100 Hz intact.
fs = 1000
t = np.arange(12 * fs) / fs
raw = (1 + .5 * np.cos(2 * np.pi * 6 * t)) * np.cos(2 * np.pi * 107 * t)
proper, u, _, _ = rate_change(np.abs(raw), fs, 100, 30, 45)
premature, v, _, _ = rate_change(raw, fs, 100, 30, 45)
proper_gain, _ = tone_estimate(proper[:, 0], u, 6)
premature_gain, _ = tone_estimate(premature[:, 0], v, 6)
assert proper_gain > .2 and premature_gain < 1e-4
checks["demodulation_order"] = {"envelope_6Hz_gain": float(proper_gain),
                               "premature_raw_reduction_6Hz_gain": float(premature_gain)}

# Independent logistic reference: training standardization, class balancing, ridge.
rng = np.random.default_rng(916)
X = rng.normal(size=(35, 3))
truth = (np.arange(35) >= 10).astype(float)
X[:, 0] += truth
X = (X - X.mean(axis=0)) / X.std(axis=0)
A = np.column_stack((np.ones(35), X))
weights = np.where(truth == 1, 35 / 50, 35 / 20)
penalty = np.diag([0, 10, 10, 10])


def objective(b):
    eta = A @ b
    return np.sum(weights * (np.logaddexp(0, eta) - truth * eta)) + .5 * b @ penalty @ b


b = np.zeros(4)
for _ in range(100):
    probability = 1 / (1 + np.exp(-np.clip(A @ b, -40, 40)))
    H = A.T @ (A * (weights * np.maximum(probability * (1 - probability), 1e-8))[:, None]) + penalty + np.eye(4) * 1e-10
    gradient = A.T @ (weights * (probability - truth)) + penalty @ b
    step = np.linalg.solve(H, gradient)
    rate = 1
    while rate > 1e-6 and objective(b - rate * step) > objective(b) + 1e-12:
        rate /= 2
    b -= rate * step
    if np.linalg.norm(rate * step) < 1e-7:
        break
reference = minimize(objective, np.zeros(4), method="BFGS", tol=1e-9)
assert abs(objective(b) - reference.fun) < 1e-8
checks["balanced_ridge_logistic_objective_difference"] = float(abs(objective(b) - reference.fun))

# Strict threshold ties, ordinal loss, and occupancy are separate quantities.
scores = np.array([.1, .2, .2, .4, .9])
cut = np.sort(scores)[::-1][int(np.floor(.2 * len(scores)))]
assert cut == .4 and np.mean(scores > cut) == .2
assert int(np.floor(.1 * 8)) == 0
checks["training_control_threshold_ties"] = True
perfect = np.eye(3)
cdf = perfect.cumsum(axis=1)[:, :2]
target = np.array([[1, 1], [0, 1], [0, 0]])
assert np.mean((cdf - target) ** 2) == 0
checks["ranked_probability_score_perfect_zero"] = True
t = np.arange(3000) / 100
g = (t >= 12) & (t < 12.5)
episode = 8 * g * np.cos(2 * np.pi * 6 * t)
assert abs(np.mean(episode ** 2) - g.mean() * np.mean(episode[g] ** 2)) < 1e-12
base = .7 * np.cos(2 * np.pi * 6 * t)
power_delta = np.mean((base + episode) ** 2) - np.mean(base ** 2)
cross_term = 2 * np.mean(base * episode)
assert abs(power_delta - np.mean(episode ** 2) - cross_term) < 1e-12
assert abs(power_delta - np.mean(episode ** 2)) > .01
checks["power_cross_term_requires_waveform_difference"] = True

# Catalog honesty and all explicit package calls/function handles resolve.
catalog = json.loads((ROOT / "methods_manifest.json").read_text())
assert len(catalog) == 112
assert sum(row["readiness"] == "implemented" for row in catalog) == 58
calls = set()
for file in ROOT.rglob("*.m"):
    source = file.read_text()
    calls.update(re.findall(r"(?:@|\b)(tbench(?:\.[A-Za-z]\w*)+)\s*(?=\(|[,;}\s])", source))
missing = []
for call in calls:
    pieces = call.split(".")
    path = ROOT.joinpath(*("+" + p for p in pieces[:-1]), pieces[-1] + ".m")
    if not path.exists():
        missing.append(call)
assert not missing, missing
checks["explicit_package_symbols_resolve"] = len(calls)
checks["implemented_catalog_entries"] = 58
checks["catalog_entries_rejected_until_supported"] = 54
checks["MATLAB_files"] = len(list(ROOT.rglob("*.m")))

(ROOT / "validation" / "study_numerical_checks.json").write_text(json.dumps(checks, indent=2) + "\n")
print(json.dumps(checks, indent=2))
