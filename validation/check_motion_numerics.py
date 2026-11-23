"""Software-4 independent numerical cross-checks and explicit source audits.

This executes Python/NumPy/SciPy reference calculations, NOT MATLAB source.
It cannot certify MATLAB runtime behavior or patient performance.
"""
from pathlib import Path
import hashlib
import json
import re
import numpy as np
from scipy.signal import fftconvolve, windows, periodogram, hilbert

ROOT = Path(__file__).resolve().parents[1]
CHECKS = {"MATLAB_native_tests_executed_here": False,
          "patient_data_analyzed_here": False,
          "scope": "independent_equations_plus_source_wiring_not_MATLAB_execution"}
RNG = np.random.default_rng(1987)


def aa_filter(fs, pass_hz, stop_hz):
    order = 2*int(np.ceil(max(8, (90-8)/(2.285*2*np.pi*(stop_hz-pass_hz)/fs))/2))
    beta = .1102*(90-8.7)
    for _ in range(12):
        n = np.arange(-order//2, order//2+1)
        h = (pass_hz+stop_hz)/fs*np.sinc((pass_hz+stop_hz)/fs*n)*np.kaiser(order+1, beta)
        h /= h.sum()
        nf = 2**int(np.ceil(np.log2(max(16384, 16*len(h)))))
        f = np.fft.rfftfreq(nf, 1/fs)
        H = np.abs(np.fft.rfft(h, nf))
        pass_error = np.max(np.abs(20*np.log10(H[f <= pass_hz])))
        attenuation = -20*np.log10(H[f >= stop_hz].max())
        if pass_error <= .05 and attenuation >= 90:
            return h, pass_error, attenuation
        order = 2*int(np.ceil(1.15*order/2))
    raise AssertionError("Anti-alias specification not met")


def fir(x, h):
    if x.ndim == 1:
        x = x[:, None]
    good = np.isfinite(x)
    y = fftconvolve(np.where(good, x, 0), h[:, None], mode="same", axes=0)
    support = fftconvolve(good.astype(float), np.ones((len(h), 1)), mode="same", axes=0)
    y[support < len(h)-.01] = np.nan
    return y


def sample(t, x, query, half=20):
    if x.ndim == 1:
        x = x[:, None]
    p = (query-t[0])/np.median(np.diff(t))
    if np.all(np.abs(p-np.round(p)) < 1e-7):
        ix = np.round(p).astype(int)
        out = np.full((len(query), x.shape[1]), np.nan)
        good = (ix >= 0) & (ix < len(t))
        out[good] = x[ix[good]]
        return out
    ix = np.floor(p).astype(int)[:, None]+np.arange(-half+1, half+1)
    d = p[:, None]-ix
    w = np.sinc(d)*np.i0(8.6*np.sqrt(np.maximum(0, 1-(d/half)**2)))/np.i0(8.6)
    w[np.abs(d) > half] = 0
    w /= w.sum(axis=1, keepdims=True)
    inside = np.all((ix >= 0) & (ix < len(t)), axis=1)
    v = x[np.clip(ix, 0, len(t)-1)]
    good = inside[:, None] & np.all(np.isfinite(v) | (np.abs(w[:, :, None]) < 1e-12), axis=1)
    out = np.sum(np.where(np.isfinite(v), v, 0)*w[:, :, None], axis=1)
    out[~good] = np.nan
    return out


def reduce_motion(x, fs, out_fs):
    pass_hz = max(1.1*9, min(30, .36*out_fs))
    h, ripple, attenuation = aa_filter(fs, pass_hz, .45*out_fs)
    z = fir(x, h)
    t = np.arange(len(x))/fs
    q = np.arange(int(np.floor(t[-1]*out_fs))+1)/out_fs
    return sample(t, z, q), q, ripple, attenuation


# Native rates are not assumed to be integer multiples of 100 Hz.
tones = []
for fs in [33.3, 60., 125., 333.]:
    out_fs = min(fs, 100.)
    t = np.arange(int(round(fs*16)))/fs
    for f in [3., 6., 9.]:
        y, u, ripple, attenuation = reduce_motion(np.cos(2*np.pi*f*t), fs, out_fs)
        good = np.isfinite(y[:, 0])
        A = np.c_[np.ones(len(u)), np.cos(2*np.pi*f*u), np.sin(2*np.pi*f*u)]
        b = np.linalg.lstsq(A[good], y[good, 0], rcond=None)[0]
        gain = np.linalg.norm(b[1:])
        phase = np.arctan2(-b[2], b[1])
        assert abs(gain-1) < .006 and abs(phase) < .002
        tones.append(dict(source_Hz=fs, output_Hz=out_fs, frequency_Hz=f,
                          gain=float(gain), phase_rad=float(phase),
                          pass_error_dB=float(ripple), stop_attenuation_dB=float(attenuation)))
CHECKS["motion_fractional_rate_phase_and_amplitude"] = tones
t = np.arange(16000)/1000
x = np.cos(2*np.pi*6*t)
x[(t > 5) & (t < 5.1)] = np.nan
y, u, _, _ = reduce_motion(x, 1000, 100)
assert not np.isfinite(y[np.abs(u-5.05) < .1]).any()
CHECKS["motion_gaps_never_interpolated_across"] = True

# Hardware clock fitting is independent of EMG waveform values or labels.
a = np.array([500., 515., 530.])
b = 42+1.0002*(a-500)
center = a.mean()
fit = np.linalg.lstsq(np.c_[np.ones(3), a-center], b, rcond=None)[0]
offset = fit[0]-fit[1]*center
assert np.max(np.abs(offset+fit[1]*a-b)) < 1e-10
assert abs(360*9*.006-19.44) < 1e-12
CHECKS["hardware_clock_affine"] = dict(offset_seconds=float(offset), rate_ratio=float(fit[1]),
                                      six_ms_phase_error_at_9Hz=19.44,
                                      default_phase_gate_degrees=15)


def vector_psd(y, fs):
    n = min(len(y), round(4*fs))
    W = windows.dpss(n, 1.5, 2, norm=2).T
    nf = 2**int(np.ceil(np.log2(max(n, fs/.125))))
    P = np.zeros(nf//2+1)
    count = 0
    for i in range(0, len(y)-n+1, max(1, n//2)):
        z = y[i:i+n]-np.mean(y[i:i+n], axis=0)
        z = np.fft.rfft(z[:, :, None]*W[:, None, :], nf, axis=0)/np.sqrt(fs)
        power = np.mean(np.sum(np.abs(z)**2, axis=1), axis=1)
        power[1:-1] *= 2
        P += power
        count += 1
    return np.fft.rfftfreq(nf, 1/fs), P/count


fs = 100
t = np.arange(3000)/fs
y = np.c_[np.cos(2*np.pi*6*t), .2*np.sin(2*np.pi*6*t), np.zeros(len(t))]
rotation, _ = np.linalg.qr(RNG.standard_normal((3, 3)))
f, p = vector_psd(y, fs)
_, rotated = vector_psd(y@rotation, fs)
error = np.max(np.abs(p-rotated))
assert error < 1e-11
band = (f >= 2.5) & (f <= 9.5)
power = np.trapezoid(p[band], f[band])
assert abs(power-.52) < .01
_, rectified = periodogram(np.abs(np.cos(2*np.pi*6*t)), fs)
freq = np.fft.rfftfreq(len(t), 1/fs)
peak = freq[1:][np.argmax(rectified[1:])]
assert peak == 12
CHECKS["signed_vector_motion_power"] = dict(rotation_invariance_max_error=float(error),
                                         estimated_guarded_band_power=float(power),
                                         expected_power=.52,
                                         rectified_6Hz_coordinate_peak_Hz=float(peak))


def harmonic_hac(y, time, frequency, fs):
    u = time-time.mean()
    H = min(3, int(np.floor(.45*fs/frequency)))
    X = np.c_[np.ones(len(u)), u]
    for h in range(1, H+1):
        X = np.c_[X, np.cos(2*np.pi*h*frequency*u), np.sin(2*np.pi*h*frequency*u)]
    inverse = np.linalg.pinv(X)
    B = inverse@y
    residual = y-X@B
    D, n = y.shape[1], len(u)
    influence = np.hstack([inverse[2:4].T*residual[:, [d]] for d in range(D)])
    C = influence.T@influence
    L = min(n//4, int(np.ceil(fs/3)))
    for lag in range(1, L+1):
        term = influence[lag:].T@influence[:-lag]
        C += (1-lag/(L+1))*(term+term.T)
    eig, U = np.linalg.eigh((C+C.T)/2)
    C = (U*np.maximum(0, eig))@U.T
    G = X[:, 2:4]
    A = np.kron(np.eye(D), G.T@G/n)
    beta = B[2:4].T.ravel()
    raw = beta@A@beta
    signed = raw-np.trace(A@C)
    variance = max(0, 4*beta@A@C@A@beta+2*np.trace(A@C@A@C))
    return raw, signed, np.sqrt(variance)


ys = y[:300]+.03*RNG.standard_normal((300, 3))
ys[:, 0] += .4*np.cos(2*np.pi*12*t[:300])
raw, signed, sd = harmonic_hac(ys, t[:300], 6, fs)
other = harmonic_hac(ys@rotation, t[:300], 6, fs)
assert abs(raw-.52) < .02 and abs(signed-.52) < .02 and sd > 0
assert np.max(np.abs(np.array([raw, signed, sd])-other)) < 1e-10
CHECKS["harmonic_conditional_HAC"] = dict(raw_power=float(raw), bias_corrected_power=float(signed),
                                        conditional_power_SD=float(sd),
                                        rotation_invariance=True,
                                        whole_search_confidence_guaranteed=False)

# Approximate oscillator witness: known signal vs noise, not probability calibration.
def oscillator(y, fs):
    N, D = y.shape
    f = np.arange(3, 9.01, .25)
    F = len(f)
    scale = max(np.sqrt(np.sum(np.var(y, axis=0))/D), 1e-12)
    z = y/scale
    off_var = max(1, np.mean(np.var(z, axis=0)))
    c, s = .995*np.cos(2*np.pi*f/fs), .995*np.sin(2*np.pi*f/fs)
    x, v = np.zeros((F, D)), np.zeros((F, D))
    p00, p01, p11 = np.ones(F), np.zeros(F), np.ones(F)
    w = np.r_[.9, np.full(F, .1/F)]
    frequencies, on = [], []
    for obs in z:
        xp = c[:, None]*x-s[:, None]*v
        vp = s[:, None]*x+c[:, None]*v
        a = c*c*p00-2*c*s*p01+s*s*p11+.002
        b = c*s*(p00-p11)+(c*c-s*s)*p01
        d = s*s*p00+2*c*s*p01+c*c*p11+.002
        innovation = obs-xp
        variance = a+.15
        k0, k1 = a/variance, b/variance
        x, v = xp+k0[:, None]*innovation, vp+k1[:, None]*innovation
        p00, p01, p11 = (1-k0)*a, (1-k0)*b, np.maximum(np.finfo(float).tiny, d-k1*b)
        prior = np.r_[w[0]*.995+w[1:].sum()*.005,
                      w[0]*.005/F+.995*(.99*w[1:]+.01*w[1:].sum()/F)]
        log_like = -.5*(D*np.log(variance)+np.sum(innovation**2, axis=1)/variance)
        off = -.5*(D*np.log(off_var)+np.sum(obs**2)/off_var)
        lw = np.log(np.maximum(prior, np.finfo(float).tiny))+np.r_[off, log_like]
        w = np.exp(np.maximum(-700, lw-lw.max()))
        w /= w.sum()
        frequencies.append(np.sum(f*w[1:])/max(w[1:].sum(), 1e-300))
        on.append(1-w[0])
    return np.array(frequencies), np.array(on)


observed = y[:1600]+.005*RNG.standard_normal((1600, 3))
ff, on = oscillator(observed, fs)
strong = (t[:1600] > 3) & (on >= .8)
assert strong.mean() > .7 and abs(np.median(ff[strong])-6) < .3
ff_rotated, on_rotated = oscillator(observed@rotation, fs)
assert np.max(np.abs(on-on_rotated)) < 1e-10 and np.max(np.abs(ff-ff_rotated)) < 1e-9
CHECKS["independent_oscillator_witness"] = dict(frequency_Hz=float(np.median(ff[strong])),
                                              supported_fraction=float(strong.mean()),
                                              vector_rotation_invariance=True,
                                              on_score_is_calibrated_posterior=False)

# Episode denominators do not erase algorithm abstentions or label unknown as zero.
state = np.zeros(1000, dtype=int)
state[100:300], state[600:] = 1, -1
alarm = np.zeros(1000, dtype=bool)
alarm[100:200], alarm[400:500], alarm[600:700] = True, True, True
positive, negative, unknown = state == 1, state == -1, state == 0
tp = np.sum(alarm & positive)
fn = np.sum(~alarm & positive)
fp = np.sum(alarm & negative)
assert tp/(tp+fn) == .5 and fp/fs == 1 and np.sum(alarm & unknown)/fs == 1
assert 60*fp/negative.sum() == 15
assert np.sum(~np.zeros(1000, bool) & positive)/positive.sum() == 1
CHECKS["episode_and_abstention_accounting"] = dict(recall=.5, missed_burden=.5,
                                                 negative_alarm_seconds_per_minute=15.,
                                                 unknown_alarm_seconds=1.,
                                                 total_abstention_recall=0.)

# Conditional spectral contribution algebra preserves destructive interference.
S = np.array([[1, .9], [.9, 1]], dtype=complex)
H = np.array([[1., -.5]])
total = np.real(np.trace(H@S@H.conj().T))
direct = np.real(np.diag(S))*np.sum(np.abs(H)**2, axis=0)
cross = 2*np.real(S[0, 1]*np.sum(H[:, 0]*H[:, 1].conj()))
excluded = (direct+cross)/total
assert abs(direct.sum()+cross-total) < 1e-12
assert np.max(direct/total) > 1 and np.min(excluded) < 0
assert np.linalg.matrix_rank(np.ones((2, 2))) == 1
CHECKS["MIMO_interference_and_rank"] = dict(model_power=float(total),
                                         direct_fractions=(direct/total).tolist(),
                                         excluded_fractions=excluded.tolist(),
                                         sum_of_contributions_is_not_constrained_to_one=True,
                                         rank_deficient_inputs_reject_attribution=True)


def fit_linear(X, y, train, ridge=1):
    A, b = X[train], y[train]
    mu, sd = A.mean(axis=0), A.std(axis=0)
    sd[sd < 1e-10] = 1
    A = np.c_[np.ones(len(A)), (A-mu)/sd]
    penalty = ridge*np.eye(A.shape[1]); penalty[0, 0] = 0
    beta = np.linalg.solve(A.T@A+penalty, A.T@b)
    return np.c_[np.ones(len(X)), (X-mu)/sd]@beta, (mu, sd, beta)


subjects = np.repeat(np.arange(12), 6)
X = np.c_[subjects/10, .1*RNG.standard_normal(len(subjects))]
y_amp = .2+subjects/8+.01*RNG.standard_normal(len(subjects))
train = subjects < 10
scale = np.median(y_amp[train & (y_amp > 0)])
pred, model = fit_linear(X, np.arcsinh(y_amp/scale), train)
modified = y_amp.copy(); modified[~train] *= 1000
after, other_model = fit_linear(X, np.arcsinh(modified/scale), train)
assert np.array_equal(pred, after) and all(np.array_equal(a, b) for a, b in zip(model, other_model))
CHECKS["held_out_target_invariance"] = dict(subjects=12, training_subjects=10,
                                         changed_test_target_factor=1000,
                                         prediction_change=float(np.max(np.abs(pred-after))))

# Compression preserves variance and PSD but can destroy when a burst occurred.
t = np.arange(3000)/100
x = np.cos(2*np.pi*6*t)*((t >= 3) & (t < 6))
shifted = np.roll(x, 1500)
assert abs(np.var(x)-np.var(shifted)) < 1e-15
assert np.max(np.abs(np.abs(np.fft.rfft(x))**2-np.abs(np.fft.rfft(shifted))**2)) < 1e-9
timing_mse = np.mean((x-shifted)**2)
assert timing_mse > .09
retained = (1.-.7)/(1.-.4)
negative = (1.-1.2)/(1.-.4)
assert abs(retained-.5) < 1e-12 and negative < 0
CHECKS["compression_capacity"] = dict(same_variance_and_global_PSD=True,
                                     disjoint_burst_waveform_MSE=float(timing_mse),
                                     example_operational_retention=retained,
                                     harmful_compression_retention=negative,
                                     total_information_in_bits_claimed=False)

# Source inspections are recorded separately from executed equations.
src = lambda name: (ROOT/name).read_text()
tune = src("+tbench/+study/tune.m")
header = re.search(r"function[^\n]*tune\(([^)]*)\)", tune).group(1).split(",")
assert len(header) == 9 and "if nargin<9" in tune and "motionViews{index}" in tune
extract = src("+tbench/+study/extract.m")
assert "taskPackets=cell(K,numel(cfg.trials))" in extract and "clear taskPackets taskShared" in extract
assert "isfield(p,'motionKey')" in extract and "p.motionKey,motionKey" in extract
alg = src("+tbench/+motion/inferenceConfig.m")
assert all("'"+key+"'" in alg for key in ["motion", "motionRuntime", "clinical", "tuning"])
for name in ["extract", "inspectSelected", "calibrate", "metrology", "temporalProbe", "controlInjection"]:
    assert "motion.inferenceConfig(cfg)" in src("+tbench/+study/"+name+".m")
loss = src("+tbench/+motion/lossRows.m")
assert "power(isfinite(power)&power<0)=0" in loss and "cfg.motion.baselineRidge" in loss
assert "intersect(trainingIds,testIds)" in loss
score = src("+tbench/+motion/scoreTrial.m")
assert "declaredTemporal" in score and "r.FNseconds=target.positiveSeconds" in score
assert "EMG_source_clock_unverified" in score
assert "m.sourceSynchronized=cfg.synchronized" in src("+tbench/+motion/task.m")
reference = src("+tbench/+motion/reference.m")
assert "s.calibrationPolicyHash,policy" in reference and "intervalScope" in reference
assert "multitaperState" in reference and "oscillatorState" in reference
policy = src("+tbench/+motion/policyHash.m")
assert "fileHash(path)" in policy and "readoutRidges" not in policy
CHECKS["source_wiring_inspected_not_executed"] = [
    "ninth tuning argument retained", "temporary per-task waveform packets",
    "motion cache identity required", "EMG hooks receive stripped config",
    "test target cannot fit readout", "common fixed-ridge tonic/MVC comparator",
    "missing physical targets stay NaN", "native dynamic abstention remains a miss",
    "negative reference requires matching calibrated implementation policy",
    "three independent reference-sensitivity rows", "conditional intervals labelled",
    "multi-muscle timing requires declared shared EMG clock"]
pdf = ROOT/"docs/tremor_motion_capture_benchmark_addendum.pdf"
assert hashlib.sha256(pdf.read_bytes()).hexdigest() == "181d7ceb7b5d6407c1d59121ce93b249aa1f8a2b05233016f8f8583b0e961412"
CHECKS["source_motion_addendum_sha256"] = hashlib.sha256(pdf.read_bytes()).hexdigest()
CHECKS["MATLAB_files"] = len(list(ROOT.rglob("*.m")))
out = ROOT/"validation/motion_numerical_checks.json"
out.write_text(json.dumps(CHECKS, indent=2)+"\n")
print(json.dumps(CHECKS, indent=2))
