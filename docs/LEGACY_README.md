# Tremor sEMG benchmark platform for MATLAB

This is a source package for comparing processing pipelines on your study. It does not contain patient data or a completed comparison of your 25 subjects. The expected complete design is 25 subjects x 13 tasks x 3 repetitions = 975 multimuscle trials / 14,625 channel records.

The numerical defaults are starting settings, not a validated clinical detector. In particular, the background ratio of 4 is an OS-CFAR-style heuristic; it is **not** a claimed false-alarm probability. Calibrate it before interpreting detection yield.

## Quick start

Unzip the folder and add its **root** to the MATLAB path. Do not add `+tbench` itself.

```matlab
addpath('C:\work\tremor_benchmark_matlab');
cfg = tbench.defaultConfig();
subjects = tbench.subjectsFromFolder('D:\study\subjects', '*.mat');

% Optional: if MVC is not a 15-element EMG-amplitude vector in each MAT file:
for s = 1:numel(subjects)
    subjects(s).mvc = MVC_by_subject(s,:); % same amplitude units and column order
end

% Replace the placeholders with your actual acquisition column order:
% cfg.muscleNames = {'...', ...}; % exactly 15 names
cfg.synchronized = true; % only if the 15 channels share a common clock
cfg.outputDir = 'D:\study\benchmark_results';

preview = tbench.inventory(subjects,cfg);
study = tbench.run(subjects,cfg);
```

`run_comparison.m` is an editable entry point. `examples/in_memory_manifest.m` adapts an existing cell array of 25 subject structures. For one subject already in the workspace:

```matlab
subjects = struct('id','S01','EMGData',EMGdata,'mvc',MVCvalues);
cfg = tbench.defaultConfig();
cfg.expectedSubjects = 1;
study = tbench.run(subjects,cfg);
```

## Exact input contract

Each subject has a task structure:

```matlab
EMGData.P4.trial2   % N x 16 real numeric matrix
```

* Columns 1:15: the 15 muscle channels, in one fixed documented order.
* Column 16: monotonically increasing time. Seconds by default; set `cfg.timeScale=0.001` for milliseconds.
* Postural tasks: P1, P2, P3, P4, P5, P6, P7.
* Kinetic tasks: K1, K2, K3, K4, K6, K7. K5 is not expected.
* Trial fields: trial1, trial2, trial3.
* `EMGData` and `EMGdata` root variables are accepted. A file with both requires `subjects(s).variable` to disambiguate.

One MAT file per subject is supported by the folder adapter. For a combined MAT file, explicitly extract its 25 subject structures and use the in-memory adapter; the code does not guess an unseen subject wrapper.

Each manifest entry requires a unique `id` and either `file` or `EMGData`/`EMGdata`. Optional `variable` selects the root MAT variable; optional `mvc` overrides a vector named `MVC` in that file. Raw MVC traces are intentionally not interpreted without an extraction rule. An absent MVC produces NaN normalized powers; absolute comparisons can still run. A force-only MVC is not a voltage normalization reference.

Sampling frequency is inferred from column 16 and optionally checked against `cfg.nominalFs`. Small integer sample gaps are expanded as NaNs and retained on the original clock. Nonuniform jitter beyond `cfg.timeTolerance` is rejected. No data are concatenated across gaps. A shared time column does not establish the hardware's synchronization accuracy; you must declare that separately.

The default raw band is 20–450 Hz and therefore requires Fs above 900 Hz. Set `cfg.rawBand` to your acquisition's justified band when needed. Incompatible bands are rejected instead of silently altered. Changing acquisition/filter bandwidth changes the observable EMG quantity and needs fresh calibration.

The code never guesses muscle names, amplitude units, amplifier rails, or clinical severity. Supply known rails with `cfg.knownRails` if available. Only known rails/nonfinite values and fully flat channels are masked; high amplitude alone is not automatically rejected.

## Requirements and validation status

The implementation targets MATLAB R2019b or newer and uses base MATLAB functions only. It includes its own FIR, Welch PSD, and DPSS construction; Signal Processing, Statistics, and Parallel Computing Toolboxes are not required. Availability under a native MATLAB runtime must be confirmed with the included self-tests.

MATLAB/Octave were unavailable in the authoring session. Package-call/structural checks and independent Python numerical checks are supplied in `VALIDATION.md`; these are **not** a claim that MATLAB integration tests ran. Before processing patient data:

```matlab
addpath('tests');
report = run_selftests();
run('examples/run_demo.m');
```

The demo is explicitly synthetic (`DEMO01`, etc.). Do not treat its output as a patient-method ranking.

## Processing methods included

| Method | Implementation | Outputs |
|---|---|---|
| `welch_cfar` | Fixed raw FIR band, rectify, low-pass; mean Welch, ordered-statistic background, signed-excess peak-area integration | Static frequency/power; no invented events |
| `hann_tracker` | Local Hann spectra, explicit off-state dynamic programming, harmonic refit in the reference envelope | Frequency track, episodes, fundamental and harmonic-inclusive reconstruction |
| `dpss_tracker` | Same tracker/refit using cached DPSS tapers | Same dynamic contract |
| `switching_harmonic_kalman` | Robust harmonic Kalman filter with two-regime probability recursion, conditional on a DPSS frequency proposal | Dynamic reconstruction and an **uncalibrated** on-state probability |

The Kalman starter is an approximate conditional-frequency model. It is not the full unknown-frequency Bayesian smoother discussed in the research review; it does not provide validated frequency credible intervals. The local harmonic amplitudes and Kalman mean reconstruction can be biased at low SNR. The calibration panel measures those failures rather than claiming unbiased power.

The reference preprocessing is a symmetric windowed-sinc FIR, applied once as a centered noncausal convolution. Its power response is `abs(H).^2`, **not** the `abs(H).^4` of a forward-backward filter. This differs from a Butterworth implementation. If you replace it, keep the reporting domain fixed across candidates and rerun transfer/calibration checks.

## Add your own method

Copy and rename `examples/custom_method_template.m`, put it on the path, then register a function:

```matlab
entry = cfg.methods(1); % copies the registry field layout
entry.name = 'my_method';
entry.run = @my_method;
entry.params = struct('myParameter',1);
cfg.methods(end+1) = entry;
```

Signature: `out = my_method(c,cfg,params)`. Start with `out=tbench.emptyOutput(c)`. The plugin receives **only one muscle**, including its amplitude envelope, detrended envelope, energy observer, validity masks, time, Fs, and MVC reference. It also receives `raw`, `rawTime`, `rawFs`, and `rawValid` for alternative demodulation or carrier-domain methods. It does not receive the other 14 channels; cross-muscle target leakage is therefore avoided by interface design. This is a single-muscle processing comparison. A future multichannel method needs a separate leave-target-group-out adapter; do not secretly pull other channels from globals.

Static methods supply scalar `frequency`, `powerMean`, `detected`, and optional spectrum. Dynamic methods additionally supply N-element `signal`, `fundamental`, `frequencyTrack`, `active`, `valid`, and capabilities. `onProbability` is optional. `tbench.validateOutput` checks the contract. Use zero reconstructed rhythm in valid off periods, NaN across invalid data, and no confident frequency in unsupported gaps. Power is in squared input-amplitude units, not PSD density or arbitrary transform magnitude.

Map a raw-domain method back onto `c.t` and the fixed reference-envelope amplitude domain before returning power/waveforms. For example, use its frequency/activity proposal with `tbench.fitHarmonics(c.y,c.t,frequencyTrack,active,c.dynamicValid,cfg)`. This measures the same envelope-modulation quantity across methods; it does not claim that neural drive or mechanical tremor power has been measured. Keep the independent benchmark observers fixed while changing the candidate's internal processing.

## Benchmark mapping and limits

| Protocol group | Executable implementation |
|---|---|
| B1: held-out modulation prediction | `predictive.m`: target-muscle energy observations, training-only source selection, 3 blocked target folds, common slow baseline and gamma working score. Source methods see only their own channel. |
| B1 temporal extension | `temporalProbe.m`, opt-in: raw masking before the complete pipeline, followed by an explicitly named **static-frequency readout adapter**. It tests masked frequency information, not native episode prediction. |
| B2: physiological corroboration | `frequencySupport.m`: other-muscle fixed-envelope/energy witnesses; `network.m`: 105 pairs and complex relationships at reference frequencies. Global witnesses are not local ground truth. |
| B3: calibration | `calibration.m`: known Gaussian AM generator, frequency/depth/intermittency/noise/drift/transient challenges; artifact-free expected reporting-envelope fundamental is the target. |
| B4: metrology | `metrology.m`: gain and polarity through raw preprocessing and each estimator. Unit/sampling/DPSS/Parseval checks also appear in self-tests. |
| B5: repeatability | `repeatability.m`: detection discordance, frequency-pair coverage, conditional positive log-power disagreement, ICC(A,1), and MVC-referenced 15-muscle compositions. |
| B6: residual / compression | `residualAudit.m`: band power and variance cross-term accounting; `predictive.m` and `aggregate.m`: target-specific prediction utility retained by each feature level. |

The primary score is held-out **cross-muscle modulation predictability**, not sensitivity/specificity or mechanical contribution. Gamma scoring is a working marginal score with training-estimated shape; temporal samples are not independent subjects. The default null panel is small and cannot establish a precise rare false-alarm probability. Increase cases and preregister challenges for a formal study.

Ordinary coherence is not an optimization objective: separate linear filters preserve population MSC at nonzero transfer. Static estimators do not get artificial coherence from fitted sinusoids; waveform-dependent graph metrics are NaN for them.

`activePower90` is the 90th percentile of **instantaneous squared reconstructed fundamental samples during active time**, not a smoothed local-amplitude quantile. Retain the waveform/events to compute another definition. `powerActive` and `powerMean` are actual reconstructed mean squares; `powerMean = occupancy * powerActive` holds for dynamic outputs with zero off-state signal, on the same valid set.

The `events` compression readout currently keeps episode start/stop, mean frequency, and power but resets its phase per episode. It therefore tests a genuinely phase-discarding event summary. `rich` retains the reconstructed fundamental and its analytic quadrature. `static` uses a frequency-conditioned sine/cosine; `yield` contains only trial-level constants, absorbed by the baseline. The resulting retention score concerns **within-record cross-muscle temporal modulation**, not all uses of power summaries or a fraction of all physiological information. Unsupported event outputs remain NaN. Negative or >1 retained-utility estimates are not clipped; tiny/nonpositive rich gains have no ratio.

## Study splitting, failures, and tuning

The registry defaults are prespecified; no cohort-learned fitting occurs by default. `foldPlan.csv` nevertheless keeps the 25 subjects grouped into five outer folds. An optional callback:

```matlab
cfg.parameterSelector = @myTrainingOnlySelector;
% selected = myTrainingOnlySelector(trainSubjectManifest,innerFoldPlan,cfg,methods)
```

receives training subject entries and an inner subject-fold plan only. It must retain method names/order and return tuned parameter sets without consulting outer-test outcomes. Fold-specific sets are retained in `study.foldParameters`. This is the extension point for nested parameter search, not an automatic claim that defaults have been optimally tuned.

Scores are summarized within subjects and given equal postural/kinetic weight. Task counts and coverage are recorded. Paired contrasts bootstrap subjects, not channels/windows. Missing trials create explicit NaN feature rows and inventory entries. Failed source methods receive baseline predictive performance with coverage zero; absent native capabilities remain unavailable. No single universal winner is fabricated.

## Main outputs

Every run creates a new timestamped subfolder, retaining earlier runs:

| File | Contents |
|---|---|
| `methodPanel.csv` | Descriptive method-level scorecard with counts and calibration axes |
| `inventory.csv`, `features.csv` | Planned records, errors, raw/normalized power, occupancy, coverage |
| `events.csv` | Subject/task/trial/muscle/method event manifest |
| `prediction.csv` | B1 folds, source channels, gamma-score gain, MSE skill, coverage, compression level |
| `temporalPrediction.csv` | Opt-in raw-masked static-adapter probe |
| `frequencySupport.csv`, `network.csv` | Witness and intermuscular relationships |
| `residualAudit.csv` | Residual band power and covariance identity checks |
| `calibration.csv`, `metrology.csv` | Known-generator and gain/polarity checks |
| `repeatability.csv`, `ICC.csv`, `topography.csv` | Repetition agreement and muscle distributions |
| `subjectScores.csv`, `pairedContrasts.csv`, `retention.csv` | Subject-level comparisons and feature-loss results |
| `foldPlan.csv`, `configuration.json`, `study.mat` | Reproducible splits, settings, parameters, and full tables |
| `signals/*.mat` | Common envelope/energy/time and individual method outputs, when enabled |

CSV `NaN` means unavailable/undefined, not a zero measurement. Empty benchmark tables state `no_rows_or_unsupported`. MATLAB tables are optional: `T=struct2table(study.features);`. Inspect a retained record with `tbench.plotRecord(file,'dpss_tracker',1)`.

At full cohort size, retained dynamic arrays may occupy gigabytes; summaries alone are much smaller. `cfg.saveSignals=false` disables them, but then you cannot reconstruct the entire information-loss audit from summary CSVs. For a quick pilot, use a few subjects/tasks before the complete 25-subject run. B1 and all-pair spectra dominate runtime. `cfg.bench.temporal=true` adds substantial rerun cost.

## Recommended first scientific run

1. Confirm time units, amplitude units, synchronized acquisition, muscle order, and MVC provenance.
2. Run self-tests and the synthetic demo in your MATLAB installation.
3. Inventory all 975 planned multimuscle trials; adjudicate missing/invalid records.
4. Pilot two subjects with all tasks. Inspect waveforms and calibration failures.
5. Freeze bands, method settings, source groups, and acceptable power/burst error tolerances.
6. Run all 25 subjects, keeping the benchmark axes and upper-tail/intermittent cases separate.

No clinical severe-subject labels, natural no-tremor ground truth, or kinematic tremor power are inferred from the supplied structure. Independent clinical/mechanical validation remains a separate question.
