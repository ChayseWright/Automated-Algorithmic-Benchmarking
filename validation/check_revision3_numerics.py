"""Independent revision-3 numerical checks; DOES NOT execute MATLAB source."""
from pathlib import Path
import itertools
import json
import re
import numpy as np
from scipy.linalg import hankel
from scipy.signal import fftconvolve

ROOT = Path(__file__).resolve().parents[1]
checks = {"MATLAB_native_tests_executed_here": False,
          "patient_data_analyzed_here": False}


def runs(mask, minimum=1):
    edge = np.diff(np.r_[False, mask, False].astype(int))
    return [(a, b) for a, b in zip(np.flatnonzero(edge == 1),
                                  np.flatnonzero(edge == -1)) if b-a >= minimum]


fs, length, rank = 100, 64, 2
t = np.arange(3000)/fs
rng = np.random.default_rng(77)
y = .03*rng.standard_normal(len(t)) + ((t >= 20) & (t < 28))*np.cos(2*np.pi*6*t)
valid = np.ones(len(t), dtype=bool)
valid[1000:1020] = False
C00 = np.zeros((length, length))
C10 = C00.copy()
all0, all1 = [], []
covered = np.zeros(len(t), dtype=bool)
for a, b in runs(valid, 40):
    v = y[a:b]-np.mean(y[a:b])
    H0 = hankel(v[:length], v[length-1:-1])
    H1 = hankel(v[1:length+1], v[length:])
    C00 += H0@H0.T
    C10 += H1@H0.T
    all0.append(H0)
    all1.append(H1)
    covered[a:b] = True
H0 = np.hstack(all0)
H1 = np.hstack(all1)
error = max(np.max(np.abs(C00-H0@H0.T)), np.max(np.abs(C10-H1@H0.T)))
assert error < 1e-10 and np.array_equal(covered, valid)
C00 /= H0.shape[1]
C10 /= H0.shape[1]
eigenvalues, U = np.linalg.eigh(C00)
U = U[:, ::-1]
grid = np.arange(3, 9.0001, .025)
steering = np.exp(2j*np.pi/fs*np.arange(length)[:, None]*grid)
stat = 1/np.maximum(np.sum(np.abs(U[:, rank:].T@steering)**2, axis=0), np.finfo(float).eps)
music = grid[np.argmax(stat)]
Us = U[:, :rank]
esprit = np.angle(np.linalg.eigvals(np.linalg.lstsq(Us[:-1], Us[1:], rcond=None)[0]))*fs/(2*np.pi)
denom = Us.T@C00@Us
pencil = np.angle(np.linalg.eigvals((Us.T@C10@Us)@np.linalg.inv(
    denom + np.eye(rank)*max(np.trace(denom), np.finfo(float).tiny)*1e-10)))*fs/(2*np.pi)
esprit = esprit[(esprit >= 3) & (esprit <= 9)][0]
pencil = pencil[(pencil >= 3) & (pencil <= 9)][0]
assert max(abs(music-6), abs(esprit-6), abs(pencil-6)) < .1
checks["full_record_subspace_late_burst"] = {
    "burst_start_seconds": 20, "burst_stop_seconds": 28,
    "valid_samples_used": int(covered.sum()), "inference_coverage_fraction": 1.,
    "pooled_covariance_equivalence_max_error": float(error),
    "MUSIC_frequency_Hz": float(music), "ESPRIT_frequency_Hz": float(esprit),
    "matrix_pencil_frequency_Hz": float(pencil)}

# Every supported sample, including tails and segments shorter than a block,
# must enter at least one atomic block. Gaps are never joined.
block_cases = 0
for total, block in itertools.product([40, 41, 160, 300, 2999], [64, 128, 192, 256]):
    mask = np.ones(total, dtype=bool)
    if total > 300:
        mask[120:140] = False
    used = np.zeros(total, dtype=bool)
    for a, b in runs(mask, 40):
        n = min(block, b-a)
        starts = sorted(set(list(range(a, b-n+1, block//2))+[b-n]))
        for first in starts:
            assert mask[first:first+n].all()
            used[first:first+n] = True
    expected = np.zeros(total, dtype=bool)
    for a, b in runs(mask, 40):
        expected[a:b] = True
    assert np.array_equal(used, expected)
    block_cases += 1
checks["atomic_full_coverage_block_tail_cases"] = block_cases


def aa_filter(fs_in, pass_hz, stop_hz):
    beta = .1102*(90-8.7)
    order = 2*int(np.ceil(max(8, (90-8)/(2.285*2*np.pi*(stop_hz-pass_hz)/fs_in))/2))
    for _ in range(12):
        n = np.arange(-order//2, order//2+1)
        h = 2*(pass_hz+stop_hz)/2/fs_in*np.sinc((pass_hz+stop_hz)/fs_in*n)*np.kaiser(order+1, beta)
        h /= h.sum()
        nfft = 2**int(np.ceil(np.log2(max(16384, 16*len(h)))))
        f = np.fft.rfftfreq(nfft, 1/fs_in)
        H = np.abs(np.fft.rfft(h, nfft))
        if np.max(np.abs(20*np.log10(H[f <= pass_hz]))) <= .05 and -20*np.log10(H[f >= stop_hz].max()) >= 90:
            return h
        order = 2*int(np.ceil(1.15*order/2))
    raise AssertionError("Unmet anti-alias specification")


def warp(v, factor):
    if factor <= 1:
        target = fs*factor
        h = aa_filter(fs, .4*target, .48*target)
        w = fftconvolve(v, h, mode="same")
        support = fftconvolve(np.ones(len(v)), np.ones(len(h)), mode="same")
        w[support < len(h)-.01] = np.nan
        v = w
    position = np.arange(int(np.floor((len(v)-1)*factor))+1)/factor
    indices = np.floor(position).astype(int)[:, None] + np.arange(-19, 21)
    distance = position[:, None]-indices
    w = np.sinc(distance)*np.i0(8.6*np.sqrt(np.maximum(0, 1-(distance/20)**2)))/np.i0(8.6)
    w[np.abs(distance) > 20] = 0
    w /= w.sum(axis=1, keepdims=True)
    inside = np.all((indices >= 0) & (indices < len(v)), axis=1)
    values = v[np.clip(indices, 0, len(v)-1)]
    good = inside & np.all(np.isfinite(values) | (np.abs(w) < 1e-12), axis=1)
    result = np.sum(np.where(np.isfinite(values), values, 0)*w, axis=1)
    result[~good] = np.nan
    return result


warp_rows = []
for factor in [1/1.9, 1/1.3, 1.3, 1.9]:
    tt = np.arange(2000)/fs
    v = warp(np.cos(2*np.pi*6*tt), factor)
    u = np.arange(len(v))/fs
    f = 6/factor
    X = np.column_stack([np.ones(len(u)), np.cos(2*np.pi*f*u), np.sin(2*np.pi*f*u)])
    good = np.isfinite(v)
    b = np.linalg.lstsq(X[good], v[good], rcond=None)[0]
    gain = np.linalg.norm(b[1:])
    assert abs(gain-1) < .005
    warp_rows.append({"warp_factor": factor, "expected_frequency_Hz": f, "gain": float(gain)})
checks["IRASA_guarded_warp_tones"] = warp_rows

# Check theory repairs and deterministic finite-space coverage separately.
for nw, k in itertools.product([1.5, 2, 2.5, 3], [1, 2, 3, 4]):
    assert k <= 2*max(nw, (k+1)/2)-1
checks["DPSS_admissible_parameter_repairs"] = 16
spec = {"highpass": [20, 35, 40], "window": [1, 1.5, 2, 3, 4],
        "threshold": [2, 3, 4, 6, 8], "minimum_cycles": [3, 4, 5]}
base = dict(highpass=20, window=2, threshold=4, minimum_cycles=3)
anchors = {tuple(base.values())}
for key, choices in spec.items():
    for value in choices:
        p = dict(base)
        p[key] = value
        p["window"] = max(p["window"], p["minimum_cycles"]/3)
        anchors.add(tuple(p.values()))
for j, choices in enumerate(spec.values()):
    assert set(choices).issubset({p[j] for p in anchors})
checks["coverage_anchor_values_represented"] = True

# Wiring audits are explicitly source inspection, not dynamic tests.
opt = (ROOT/"+tbench/+study/optimize.m").read_text()
assert "allowed=initialIndices" in opt and "find(training),trialIndices" in opt
assert "selectionScope" in opt and "all_training_deployment" in opt
extract = (ROOT/"+tbench/+study/extract.m").read_text()
assert "candidateBatchSize" in extract and "preprocessOwner" in extract
cross = (ROOT/"+tbench/+study/validatePair.m").read_text()
assert "requires %s transform coefficients" in cross and "prepared DPSS replicates" in cross
checks["source_wiring_inspected_not_executed"] = [
    "fold-local adaptive pool", "training-only proposal extraction", "explicit deployment scope",
    "bounded candidate batches", "preprocessing-owner cache identity", "transform compatibility"]
checks["MATLAB_files"] = len(list(ROOT.rglob("*.m")))
path = ROOT/"validation/revision3_numerical_checks.json"
path.write_text(json.dumps(checks, indent=2)+"\n")
print(json.dumps(checks, indent=2))
