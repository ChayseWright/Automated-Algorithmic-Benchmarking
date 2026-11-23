# Signal model, reduction and power definitions

## What information is being preserved

Surface EMG is a superposition of motor-unit action potentials with amplitude cancellation, recruitment and discharge modulation. Calling the complete signal a conventional narrowband FM carrier is generally inadequate. Low-frequency rhythmic modulation can be observed in a rectified envelope or carrier-energy statistic, but these are nonlinear proxies for neural drive.

For a textbook signal

x(t) = A(t) cos(2*pi*fc*t + beta*sin(2*pi*fT*t)),

a pure phase/frequency modulation with constant A can leave the amplitude envelope nearly constant. Therefore preserving an envelope does not preserve every possible FM observable. A separate carrier branch is retained for methods that need it. The FM audit checks preservation of a known waveform under carrier reduction; it does not establish that a particular algorithm detects ET from pure FM.

The physiological limitation is supported by Dideriksen and Farina (2019), doi:10.1371/journal.pcbi.1006985, [PubMed](https://pubmed.ncbi.nlm.nih.gov/31050667/). See also the unchanged protocol PDF.

## Fixed reporting observer

For every trial, at its original measured acquisition rate:

1. Apply one centered symmetric Hamming-windowed 20–450 Hz FIR (0.15 s nominal order support).
2. Rectify for the reporting envelope; square for a separate energy observer.
3. Apply a measured-response anti-alias FIR.
4. Sample both observers on the same 100 Hz grid.
5. Subtract a 1 Hz slow component for oscillatory harmonic fitting; retain the undetrended envelope and positive energy separately.

The anti-alias passband is 0–30 Hz, stopband begins at 45 Hz, declared attenuation is at least 90 dB and pass error is at most 0.05 dB. The filter order is increased until its sampled response meets both limits. The 27 Hz third harmonic of 9 Hz lies in the preserved passband.

Centered convolution compensates FIR delay but is noncausal. Invalid sample support and boundaries remain invalid. This is an offline 30-second estimator, not a real-time latency specification.

The reporting observer is computed before any raw carrier reduction. A high-pass choice made for detection does not redefine the measured envelope or MVC power denominator.

IRASA intentionally stretches/compresses sample sequences to move oscillatory spectral peaks relative to a scale-free background. Its warps therefore use a separate `preprocess.irasaWarp` module: down-warps are anti-aliased, up-warps use guarded sinc interpolation, and spectra are interpreted at the original analysis rate by design. These are an analysis operation, not a claim to preserve the physical waveform unchanged.

## Two reduction branches

| Branch | Default behavior | Reason |
|---|---|---|
| Envelope/energy | Original-rate nonlinearity, anti-alias filter, then 100 Hz | Preserve target modulation and declared harmonics |
| Carrier | Integer factor floor(Fs/1200), never below one | Keep a separate raw representation through the declared 450 Hz band |
| Original MAT | Unchanged | Reproduce other observables and investigate discarded bands |

For example, 2,500 Hz gives a 1,250 Hz carrier branch and 100 Hz observers. A 1,111 Hz carrier remains at its original rate. Reduction-enabled controls the carrier branch; the envelope reporting clock remains the separately declared `analysisFs`. To study envelope-rate sensitivity, set a higher `analysisFs` while retaining the same declared envelope pass/stop edges and rerun a frozen experiment.

Carrier anti-alias passband ends at 450 Hz; the stop edge is placed within its output Nyquist margin. It is not a raw-to-100 Hz shortcut. Information outside the declared carrier band is not guaranteed to survive that branch.

An integer clock ratio uses direct sampling after FIR filtering. A fractional ratio uses normalized 40-tap windowed-sinc interpolation after filtering, preserving the time grid across channels. NaN neighborhoods remain unsupported. The interpolation tolerance and anti-alias response are separate: a measured FIR response alone is not proof of ideal fractional interpolation.

[MATLAB resample documentation](https://www.mathworks.com/help/signal/ref/resample.html) describes the standard FIR anti-alias / rate-change principle. This package supplies its own base-MATLAB guarded implementation; it does not call that toolbox function or claim bitwise equivalence to it.

## What “no information loss” can and cannot mean

The original files are not replaced. No selected waveform/evidence is converted to single precision. Repeated clocks and common observer arrays are stored once and restored by `tbench.loadSelected`. Those storage changes avoid duplication without changing values.

Filtering/resampling is not exactly invertible. Losses include out-of-band content, finite response ripple, interpolation approximation and unsupported boundaries/gaps. Abrupt burst onsets contain broadband information above 30 Hz, so onset detail is not mathematically lossless even when the tremor carrier frequency is 3–9 Hz. Calibration reports timing and power consequences.

The fixed reporting domain already excludes carrier frequencies outside 20–450 Hz. This is an explicit measurement definition, not a claim that no physiological content exists elsewhere.

The independent Python checks measured the following reference implementation values:

| Input rate | Anti-alias taps | Pass error | Stop attenuation |
|---|---:|---:|---:|
| 1,000 Hz | 441 | 0.000285 dB | 91.10 dB |
| 1,111 Hz | 489 | 0.000289 dB | 91.26 dB |
| 2,500 Hz | 1,097 | 0.000293 dB | 91.61 dB |

The pure-FM 2,500→1,250 Hz reference challenge had RMS error about 1.38e-6 on supported samples. These are Python equation checks, not native MATLAB execution or patient measurements. The MATLAB suite regenerates its actual-rate audit in `transferAudit.csv`.

## Power estimators

All primary reported values use the same fixed electrical-envelope units:

- Static methods propose a frequency and accepted-record flag. A fixed-reference Welch meter integrates **signed PSD excess over local background** within ±0.75 Hz of the proposal, then constrains the displayed power to be nonnegative. The signed value remains available. The reporting Welch window is fixed by study configuration; it is separate from a candidate's detector window.
- Dynamic methods propose frequency and active support. Local harmonic regression on the fixed reference reconstructs the fundamental and harmonics. Fundamental mean power is mean(fundamental²) over valid support. Active power, occupancy and episode power are separate values.
- The method's internal power and frequency remain available as native audit values. A wavelet/coherence magnitude is not silently relabelled as electrical power. Kalman proposals plus common refitting are not a native posterior-expected-power estimate.

Static peak-area and dynamic reconstruction are different estimators. A fixed spectral band can truncate broad burst sidebands; dynamic refits can attenuate short episodes or fit noise. Known-input power/timing calibration is therefore required beside clinical scores. A common unit alone does not make their biases identical.

For valid dynamic support, mean power = occupancy × active power when activity exists and the reconstruction is zero in off periods. If no activity exists, active power is undefined, not a physiological zero.

MVC normalization is 10,000 * power / MVC², in (%MVC)². MVC remains an amplitude reference; its extraction convention must be stated. Voltage², V²/Hz, dimensionless coherence and (%MVC)² are not interchangeable.

Do not add component powers without checking cross terms. Power(estimate after injection) minus power(estimate before injection) is not generally incremental component power. The control-counterfactual adapter uses waveform differences and a declared projection where those waveforms exist.

