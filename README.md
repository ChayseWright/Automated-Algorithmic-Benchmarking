# Essential tremor sEMG benchmark suite — software 4.0

This update adds the supplied motion-capture benchmark addendum to the existing MATLAB study with ET participants, controls and TETRAS. It preserves **preprocessing → analysis → benchmark**, keeps the EMG estimators independent of motion, and retains **lower-severity capture** as the primary selection objective under control, severe-case, fidelity and coverage constraints.

**Validation status:** independent Python equation checks and structural checks passed. MATLAB/Octave is unavailable in the authoring environment; the MATLAB files and integration tests have **not** been executed here. No participant recordings or TETRAS table were supplied. There is no estimated patient result or demonstrated winning method in this package.

## Start here

Requirements: MATLAB R2019b or later with the JVM enabled. The implementation uses base MATLAB, including its optimization and linear algebra functions; Signal Processing Toolbox is not required. Octave is not a supported/tested target.

After extracting the ZIP, enter this directory in MATLAB:

~~~matlab
nativeReport = START_HERE('validate');
demo = START_HERE('demo');     % Clearly labelled synthetic data only.
motionDemo = START_HERE('motion_demo');
% Run your study only after investigating any validation failure:
study = START_HERE('run','/absolute/subject_mat_folder', ...
                  '/absolute/clinical.csv','/absolute/results');
% Separate, exact-ID-matched canonical motion MAT files:
study = START_HERE('motion','/absolute/subject_mat_folder', ...
                  '/absolute/clinical.csv','/absolute/results','/absolute/motion_folder');
% Omit the final folder when MotionData is embedded in the subject MAT files.
~~~

Start with `examples/run_study_demo.m` to inspect the output structure using clearly labelled synthetic data. Its reduced coverage settings and disabled calibration gate are demonstration settings, not the full-study defaults.

For configuration:

~~~matlab
cfg = tbench.studyConfig();
cfg.nominalFs = 1111;             % Optional check; use the actual acquisition rate.
cfg.timeScale = 1;                % Seconds; use .001 only for millisecond timestamps.
cfg.synchronized = true;          % Only when all 15 channels share an acquisition clock.
% cfg.muscleNames = {...};        % Exactly 15 acquisition-order names.
cfg.outputDir = '/absolute/results';
cfg.cacheDir = fullfile(cfg.outputDir,'cache');
study = run_benchmark('/absolute/data','/absolute/clinical.csv',cfg);
~~~

`run_comparison.m` is an editable EMG-only script. The legacy runner remains available through `tbench.defaultConfig()`; this study uses `tbench.studyConfig()`. Software 4 retains benchmark protocol version 2 and includes both that PDF and the supplied motion addendum unchanged.

## Motion-capture extension

Read [docs/MOTION_CAPTURE.md](docs/MOTION_CAPTURE.md) for the complete input schema, clock model, calibration requirements, endpoint definitions, capacity limits and exact output names. `examples/make_motion_trial.m` gives an import template; `tbench.motion.fromPoses` adapts explicitly named 6-DOF segments. Vendor C3D/CSV/marker conventions are not guessed.

The motion extension adds M1 frequency agreement, M2 episode timing, M3 missed motion burden, M4 physical waveform/RMS prediction, M5 coherence with guarded conditional attribution, and M6 unsupported mechanical alarms. Separate compression reports quantify target-specific retained utility. Static PSD methods retain their trial-level endpoints; native waveform and timing scores remain unavailable.

Motion is off in `tbench.studyConfig()` for compatibility. `run_motion_benchmark` enables it and requires at least one imported trial. Missing or uncertain motion does not remove a participant from the EMG clinical benchmark. Reference construction consumes signed motion coordinates and clock metadata, never EMG values or diagnosis; EMG hooks receive a configuration stripped of motion and clinical labels. Two additional reference witnesses expose sensitivity of frequency/episode conclusions.

All outer predictions remain subject-held-out. A low-capacity population readout predicts physical trial RMS amplitude. A separate person-calibrated waveform readout trains on two repetitions and evaluates the third. Its scope is labelled explicitly; it does not estimate unseen-person waveform generalization. Controls can have real physiological motion: mechanical negatives require an externally calibrated amplitude limit, not control membership.

## Input contract

One MAT file per participant, with a scalar variable named `EMGData` or `EMGdata`:

~~~matlab
EMGData.P4.trial2 = A;    % A is N-by-16, real numeric.
% Columns 1:15 = the same muscles in the same order.
% Column 16 = strictly increasing time.
~~~

Defaults are P1–P7, K1–K4/K6/K7, and trial1–trial3: 39 trials and 585 muscle records per participant. The assumed 25 ET plus 10 controls therefore give 1,365 multimuscle trials and 20,475 muscle records. Actual counts and valid durations are retained; missing records produce explicit failure rows.

Sampling frequency comes from the time column. Missing samples are expanded to NaNs on the original clock. Excess timestamp jitter is rejected rather than silently warped. Inf samples and values at declared ADC rails are invalidated. A raw rate incompatible with the declared 20–450 Hz carrier band is rejected for that record; the software does not silently change the physical band.

Supply MVC as either a top-level `MVC` numeric 15-element amplitude vector in each MAT file, `subjects(s).mvc`, or all relevant `mvc_01`–`mvc_15` CSV columns. These references must use the same amplitude units as the EMG. Raw MVC traces require a separately specified plateau/extraction rule. A CSV override containing some missing references leaves those references missing; it does not fill them by guessing.

The clinical CSV joins **exactly by MAT filename stem / subject ID**, never by row ordering. See `examples/clinical_template.csv`. Required columns:

| Column | Meaning |
|---|---|
| subject | Exact participant ID / filename stem |
| group | ET or control |
| severity_score | Numeric TETRAS target; blank/NaN when absent |
| severity_kind | Recommended: exact item/subscale, percent-of-maximum or percentile definition |
| severity_group | mild/moderate/severe when using provided categories |
| phenotype | Independently documented ET, ET-plus, ET_PD, unknown, etc. |

Optional side, visit and age fields are preserved. All original CSV columns/values are retained in the participant metadata, including additional version, rater and task-mapping columns. They are not automatically turned into model covariates.

Controls may have missing TETRAS. If using provided categories:

~~~matlab
cfg.clinical.severityMode = 'provided';
~~~

Alternatively, prespecify two justified cuts in `cfg.clinical.severityCuts`. Without such definitions, the software uses **training-only relative ET tertiles**, keeping equal scores in the same category. “Mild” and “severe” metric column names then mean relative lower and higher severity, not transportable clinical grades. Percentile rank is not a physical tremor amplitude.

## Three modules

| Module | Entry | Contract |
|---|---|---|
| Data preprocessing | `pipeline.preprocess(c,cfg,pipeline)` | One muscle, matched carrier/envelope/spectrum/TF representation, masks and clock |
| Analysis | `pipeline.analyze(prepared,cfg,pipeline)` | Frequency, supported power, native capabilities, optional episodes, evidence and audit |
| Benchmark | `pipeline.benchmark(predictions)` | Subject-level control/capture/severity metrics |

`tbench.pipelineRegistry` matches inputs and analysis automatically. The six primary pipelines are:

| Analysis | Matched input |
|---|---|
| Welch/OS-CFAR comparator | Rectified reference envelope → Hann Welch spectrum |
| Hann off-state tracker | Local Hann spectra |
| DPSS off-state tracker | Local DPSS spectra |
| Switching harmonic Kalman proposal | Detrended envelope with prepared local DPSS evidence for its frequency seed |
| Morse CWT proposal/tracker | Complex Morse transform |
| Multiresolution STFT tracker | Envelope with long/short local spectra |

The saved broader prototype contained entry-point wrappers without several family backends. The corrected catalog has **58 implemented research entries** and **54 unavailable entries**. Unavailable algorithms fail preflight; they are never scored as zero-power successes. `docs/method_coverage.csv` lists every entry. “Implemented” means source exists, not native validation or equivalence to a published reference toolbox. The six-pipeline comparison is the default; the broader catalog is exploratory.

## Parameter selection

Each method has a declared, method-specific search space. Small spaces are exhausted. Larger spaces start with the baseline and main-effect anchors covering every declared parameter value, then fill a reproducible budget of 64 signal configurations. Two training-only refinement rounds propose up to 16 additional configurations each, using up to three leading configurations. Clinical ridge values are evaluated without rerunning signal extraction. A search budget too small to cover the declared values is rejected unless you explicitly enable demonstration-style truncation.

Five outer subject folds estimate the **learning procedure**. Inner subject folds select parameters, transforms/readouts and the empirical control operating threshold. All tasks, trials and muscles from a person stay together. There is one outer prediction per participant, method and feature level.

Selection maximizes lower-severity capture among candidates satisfying prespecified control, severe-capture, coverage and calibration requirements. With motion enabled, tied candidates use inner-validation RMS-prediction gain when it is evaluable for every tied candidate; remaining ties use balanced log loss and the prespecified candidate order. Outer-test motion never selects parameters. No qualifying candidate produces `constraints_unmet` and an explicitly provisional least-violation choice. See `docs/PARAMETER_SELECTION.md`.

The final all-training configuration per method is saved separately for future inference. It is the **best tested** configuration under the declared rule, not a global optimum, and does not validate itself.

Use `tbench.crossPipelines` to create explicitly labelled compatible preprocessing/analysis comparisons. Transform ownership is frozen in these cross experiments; incompatible transform/container pairs are rejected. See `examples/cross_modules.m`. A matched front end is a defensible starting point, not proof that preprocessing is globally optimal.

## Reduction without changing the reporting endpoint

A low tremor frequency is not permission to decimate raw sEMG to 18–20 Hz. The fixed reporting envelope is computed **before** carrier reduction. It uses original-rate carrier filtering, rectification, measured anti-alias filtering, and a shared 100 Hz clock. The passband extends through 30 Hz to retain three harmonics of 9 Hz. Methods requiring carrier information receive a separate conservatively reduced raw branch.

Original MAT files are never overwritten. No amplitude quantization is introduced in selected waveforms or evidence. Filtering and resampling are not mathematically lossless: their declared transfer tolerances, boundary/gap exclusions and out-of-band losses are audited. See `docs/DATA_REDUCTION.md`.

Extraction processes at most 32 candidates per batch by default, bounding retained feature structures at the cost of extra source-file reads. Shared transforms and content-validated caches reduce repeated work. MUSIC, ESPRIT and matrix pencil use all supported continuous portions of the record. Atomic-norm proposals use overlapping blocks covering the full supported record and expose their finite-iteration residuals. IRASA warps use anti-alias filtering and guarded sinc sampling.

## Find results

Each run receives an atomic run number plus timestamp/unique ID. Open its `index.html`.

| File / directory | Contents |
|---|---|
| selected_parameters.csv | Per-outer-fold choices and separate all-training choices |
| candidate_catalog.csv / tuning_scores.csv | Tested configurations, hashes, constraints and inner results |
| parameter_spaces.csv / parameter_stability.csv | Declared choices and descriptive outer-fold choice frequencies |
| optimizationTrace.csv | Training membership, refinement rounds, candidate counts and termination |
| clinicalPredictions.csv / clinicalScorecard.csv | Outer scores, decisions, severity, coverage and metrics |
| rhythmSummary.csv | Raw accepted rhythm yield, including secondary 4–8 Hz yield |
| calibrationTraining.csv / calibrationTest.csv | Development and independently seeded signal challenges |
| transferAudit.csv / reduction.csv | Anti-alias, phase/FM checks and actual acquisition-rate audits |
| clinicalRetention.csv / physiologyRetention.csv | Operational information-retention comparisons |
| phenotypeScorecard.csv / qualitySummary.csv | Phenotype, task, failure and out-of-band witness summaries |
| features.csv / events.csv | Selected per-record estimates and native episode details |
| references / selected | Shared reporting signals and selected method outputs |
| models / deployment | Outer fits versus final all-training policies |
| participants.json / run_manifest.json / study.mat | Labels/metadata, provenance and native configuration |
| explanatory | Streamed detailed physiological checks |
| checkpoints | Completed-fold checkpoints; feature caches support efficient reruns |
| motion_referenceAudit.csv / motion_coverage.csv | Reference states, quality, clock, limit calibration and missing trials |
| motion_agreement.csv / motion_agreementSummary.csv | M1–M3/M6, reference sensitivity and severity/task strata |
| motion_predictions.csv / motion_waveformPrediction.csv | Separate unseen-subject RMS and person-calibrated waveform outcomes |
| motion_retention.csv / motion_waveformSummary.csv | Compression loss, operational retention and evaluability |
| motion_coherence.csv | M5 selected-method coherence and conditional contribution eligibility |
| motion_clinicalScorecard.csv | Independent motion-only clinical comparator |
| motion/references / motion/models | Shared physical reference arrays and training-fitted RMS readouts |

Load a selected result with `tbench.loadSelected(filename)`; it resolves its shared reference and restores repeated time axes. Inspect it with `tbench.plotRecord(filename,pipelineId,muscleIndex)`.
For motion-enabled records, use `tbench.plotMotionRecord(filename,'hand_torso',[5 10 14])`, substituting your actual target name and acquisition indices.
Inspect candidate tradeoffs using `tbench.plotMotionSelection('/results/run_directory/study.mat',pipelineId,outerFold)`; it shows lower-severity capture versus inner mechanical gain and marks the eligible Pareto frontier.

Apply a frozen all-training policy using `tbench.predictSubject(modelFile,subjectSpec)`. Its selection status remains explicit. It checks source/module versions and does not retune on the new participant.

The default calibration has 75 cases: 3/6/9 Hz, depths 0/0.05/0.10/0.25/0.80, and continuous/early/middle/late/repeated bursts. Low-depth frequency recovery and late strong-burst recovery have separate audit columns. Weak synthetic modulation is a controlled signal challenge; it is not a clinical definition of mild ET.

## What remains a scientific limitation

A control rhythm is not necessarily a false physiological signal; the clinical false-positive unit is a participant. Ten controls provide coarse empirical specificity estimates. EMG envelope power is not identified mechanical tremor power or latent neural drive. Global TETRAS is not a muscle-by-muscle episode annotation. The benchmark distinguishes these targets.

Compression reports are target-specific losses under defined readouts, not a percentage of all Shannon information. Phenotype tables evaluate subsets of an inclusive fit; a restricted-cohort sensitivity analysis requires a separately frozen manifest/rerun. Age and acquisition confounding are not automatically corrected.

Read `docs/BENCHMARKS.md` for exact capacities/denominators, `docs/MODULES.md` for extension rules, and `VALIDATION.md` for executed versus unexecuted checks. The source PDF is included unchanged.
