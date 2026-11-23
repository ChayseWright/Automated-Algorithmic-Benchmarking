# Motion benchmark implementation — software 4.0

This code implements the supplied `tremor_motion_capture_benchmark_addendum.pdf` as a separate mechanical benchmark panel. The original EMG clinical benchmark remains active. There are no patient results or fitted patient parameter sets in this distribution.

## Run and configure

From the extracted project directory in MATLAB:

```matlab
addpath(pwd); addpath(fullfile(pwd,'tests'));
nativeReport = run_native_validation(true);
demo = START_HERE('motion_demo');

cfg = tbench.studyConfig();
cfg.outputDir = '/absolute/results';
cfg.cacheDir = fullfile(cfg.outputDir,'cache');
cfg.motion.primaryTargets = {'hand_torso'};  % Your canonical target name(s).
cfg.motion.clinicalTargetKind = 'translation';
cfg.synchronized = true;                   % Only after verifying the 15 EMG channels share a clock.
study = run_motion_benchmark('/absolute/EMG','/absolute/clinical.csv', ...
                             '/absolute/motion',cfg);
```

Use an empty motion-folder argument for embedded `MotionData`. An in-memory manifest may supply `subjects(s).MotionData` and `subjects(s).EMGData`; then set `cfg.motion.enabled=true`, optionally `required=true`, and call `tbench.runStudy`.

`START_HERE('motion',EMGFolder,clinicalCSV,resultsFolder,motionFolder)` is the convenience entry. EMG-only `START_HERE('run',...)` retains the previous behavior. Motion settings are in `tbench.motionConfig`; preflight rejects incompatible sampling/band/window settings. The synthetic demo reduces search budgets and disables calibration gates for demonstration only.

## Canonical subject/task/trial schema

One separate motion MAT file is named by the exact subject ID, for example `ET_01.mat`, containing `MotionData`; or include that variable in the EMG subject file. The hierarchy matches EMG exactly:

```matlab
MotionData.P4.trial2 = record;
record.time = motionTimeSeconds(:);       % Native timestamps, not a resampled guess.
record.reviewed = true;                   % Only after kinetic ambiguity review.
record.sync.motionEvents = [500;530];      % Matched acquisition events, seconds.
record.sync.emgEvents = [42;72];
record.sync.uncertaintySeconds = .001;     % Bound covering the trial and event errors.

target.name = 'hand_torso';
target.kind = 'translation';
target.values = handInThoraxCoordinates;  % N-by-3 SIGNED coordinates.
target.unit = 'mm';
target.frame = 'thorax';
target.valid = motionQualityMask(:);      % N-by-1 or N-by-D finite 0/1 flags.
target.confounded = ambiguousMovement(:); % N-by-1 finite 0/1 flags.
target.muscles = 1:15;                    % Frozen source set for the limb alarm/readout.
record.targets = target;
```

For multiple targets, use a uniform struct array with the same fields. Names must be unique character rows and consistent across the cohort. The reader rejects a named target that changes coordinate kind, unit after normalization, or frame during readout fitting. IDs are joined by filename, never by row order or approximate string matching.

| Field | Supported meaning |
|---|---|
| time, timeScale | Increasing native time; optional positive scale converts time to seconds |
| values | One to three real signed coordinates; missing values are NaN |
| kind | `translation`, `rotation`, or `joint_angle` |
| unit | Translation `m`/`mm` normalized to mm; angles `deg`/`rad` normalized to radians |
| frame | Explicit physical coordinate frame; `unspecified` is preserved, not inferred |
| valid | Optional supplied hardware/geometry quality; otherwise finite-value support only, labelled `qualitySupplied=false` |
| confounded | Ambiguous rhythmic movement; excluded from binary rhythm states |
| reviewed | Target override or trial flag; kinetic binary reference states require review by default |
| muscles | Unique acquisition indices 1–15; default all 15, without an anatomical source claim |
| noiseSD, limitRMS | Optional external noise calibration and vector RMS amplitude limit in native physical units |
| calibrationVerified, calibrationPolicyHash | Explicit calibration declaration tied to the exact reference policy/implementation |

Joint angles are unwrapped separately within each continuous valid run. Quaternion/marker matrices must first pass through a declared geometry adapter. Never rectify a coordinate or take its magnitude before spectral estimation: doing so can replace a 6-Hz line by a 12-Hz component. Vector PSD uses the sum of signed-axis auto-spectra, preserving rotation invariance within one physical unit/frame. It does not add mm² and rad².

`examples/make_motion_trial.m` is an editable template. `tbench.motion.fromPoses(time,segments,referenceName,handName,unit,jointPairs)` accepts positions and active `[w x y z]` quaternions. It propagates parent/child validity, expresses hand translation in the reference frame, and optionally returns relative rotation vectors. Those vectors are **not anatomically calibrated shoulder/elbow/wrist DOFs**; values near the rotation-log singularity at pi are masked. Segment conventions, joint calibration and vendor import remain your explicit responsibility.

## Synchronization and sampling

The affine mapping is `EMG absolute seconds = offset + rateRatio * motion absolute seconds`. At least two hardware events fit offset and rate; one event also needs a finite nonnegative `clockPPMbound`. Alternatively, declare `sync.verified=true` with explicit `offsetSeconds`, `rateRatio` and `uncertaintySeconds`. Events are always supplied in seconds, including when `timeScale` converts the record timestamps from another unit.

The software never estimates acquisition offset by maximizing EMG–motion correlation. Unverified clocks retain separate trial origins: independent trial RMS/frequency remains possible; matched timing, phase, coherence and waveform scores are unavailable. Verified timing requires uncertainty <=20 ms. Phase requires `360*9*uncertainty <=15 degrees` by default, about 4.63 ms. A physiological EMG-to-motion delay is not an acquisition-clock error. The limited quadrature readout can represent a constant phase relationship but does not identify its cause.

Multi-muscle phase/timing endpoints additionally require the existing `cfg.synchronized=true` declaration. The default is false: shared matrix columns alone are not a hardware synchronization certificate. This gate does not remove static frequency, continuous RMS or the EMG clinical benchmark. Single-source endpoints use that source's verified motion alignment.

Motion is reduced to `min(native Fs,100 Hz)`, never upsampled to invent bandwidth. Default minimum native rate is 25 Hz. Anti-alias pass/stop edges adapt to that actual rate while preserving the 3–9 Hz band; measured rejection must meet 90 dB and pass error <=0.05 dB. Integer and guarded sinc paths preserve clocks and masks. Missing native samples are expanded; excessive jitter or gap expansion is rejected. No interpolation bridges invalid support. Transfer details and available harmonic bandwidth are saved per target. Use original coordinates for any future endpoint requiring discarded bandwidth.

Verified clocks restrict the mechanical endpoint to the nominal EMG trial interval. Reference `validFraction` is valid mechanical seconds divided by nominal EMG trial seconds, rather than a percentage of only the available motion fragment. Recorded, valid and nominal durations and capture fraction are exported separately. A short clean motion fragment cannot silently masquerade as a complete 30-second reference.

EMG reduction remains the original-rate carrier-filter/demodulation pipeline with its separate preserved-carrier branch. sEMG is not assumed to be a pure FM carrier. Its envelope benchmark and a carrier-information method have distinct input contracts. A small mechanical tremor frequency does not justify direct 20-Hz raw EMG decimation.

## Independent motion reference

`tbench.motion.reference` accepts motion plus EMG clock metadata only. It removes a declared slow component below 1.5 Hz, computes signed vector spectra with long/short DPSS windows (3/1 s, 0.2-s hop, NW=1.5, two tapers), and estimates the selected fundamental by harmonic regression with conditional HAC uncertainty. The local fit includes nuisance harmonics rather than assigning them all to the fundamental.

A second motion-only damped switching oscillator bank provides a distinct frequency/on-state witness. It uses one vector scale, not per-axis whitening, preserving geometry. Its on score is an approximate switching-model score, not a calibrated probability of tremor. The default consensus reference requires both witnesses to support a frequency-consistent rhythm and the conditional lower power interval to exceed zero. Minimum episode duration and cycle count are enforced.

| State | Meaning |
|---|---|
| +1 | Model-supported mechanical rhythm under the frozen reference policy |
| 0 | Indeterminate, ambiguous, unsupported, or between the positive/negative criteria |
| -1 | Below the declared calibrated amplitude limit, under the reference model |

This is an independent kinematic reference, not a diagnosis or a per-muscle activation label. No motion can coexist with isometric EMG tremor. Controls can have genuine physiological motion. Mechanical positives are therefore evaluated separately from ET/control clinical decisions.

Negative states require positive finite `noiseSD`, positive `limitRMS`, `calibrationVerified=true`, and `calibrationPolicyHash=tbench.motion.policyHash(cfg.motion)`. Set these only after external calibration of the actual system, geometry and reference pipeline. The hash binds settings and reference source files; it does not perform or prove hardware calibration. NoiseSD is calibration provenance, not a claim that all coordinate noise is independent Gaussian. Controls alone do not certify sensor noise. Missing/stale calibration leaves negatives indeterminate; M6 can then be unavailable.

Intervals are conditional on a selected frequency and model. They do not provide a whole-search false-alarm or confidence guarantee. `motion_calibration.csv` is a synthetic software challenge, not instrument certification. It covers 3/6/9 Hz, null/low/high physical amplitudes, continuous/late bursts and occasional gaps. Review its false positives and coverage rather than interpreting the synthetic calibration flag as empirical device validation.

The consensus, multitaper-only and oscillator-only positive-state/frequency witnesses are scored separately in `motion_agreement.csv`. The oscillator alone does not certify negative amplitude-limit intervals. All retain the same signed-coordinate continuous power target; this is a reference-sensitivity analysis, not three independent ground truths. Kinetic states are indeterminate until review by default. Continuous RMS targets remain measured band-limited motion, potentially including voluntary rhythmic movement, rather than automatically being called ET tremor.

## Endpoint capacities and denominators

| Panel | Implemented statistic | Capacity and scope |
|---|---|---|
| M1 | Trial frequency error; missing estimate; dynamic overlap-conditional tracking error; harmonic flag | Trial estimate available to static methods; dynamic within-tolerance recall and track-missing fraction retain all positive time |
| M2 | Limb-alarm recall, resolved precision, IoU; matched-event onset/offset error | Native dynamic methods and reliable common clock; unknown reference intervals are not negatives |
| M3 | Fraction of supported positive local motion burden while the limb alarm is off | Pair burden weighting with unweighted mild-subject recall |
| M4a | Physical vector waveform MSE and gains over zero and a common broad-envelope readout | Person-calibrated: two repetitions fit, third evaluates, rotated three ways |
| M4b | Unseen-subject trial RMS error and transformed prediction gain over tonic effort/MVC | Subject-grouped outer testing, low-capacity common linear readout |
| M5 | Vector-output trace-weighted ordinary coherence; multiple coherence; conditional source contributions | Selected outputs only; rank, condition, replicate and model-fit eligibility gates |
| M6 | Positive limb-alarm seconds per reference-negative minute | Requires calibrated negative intervals and native temporal output |

The frozen limb alarm defaults to any native active output among the declared target muscles, with >=80% source coverage. An algorithm failure/abstention remains a miss during evaluable positive time; source coverage is reported alongside it. Static methods have unavailable native timing/waveform endpoints, not invented zero error or zero false-positive time. Frequency disagreement near 2:1 is retained as a harmonic relation rather than automatically “corrected.”

Onset errors use one-to-one event matches and are conditional on successful matching. Episode counts and recall expose missing episodes. Nominal boundary errors and excess beyond reference-window/clock uncertainty are both retained. Whole-window reference uncertainty is not silently treated as millisecond event truth.

Waveform training uses only signed method waveforms and their quadratures, never motion-selected frequency/phase as input features. The default training stride retains at least four samples per highest primary-band cycle. It has at most two columns per declared EMG source plus an intercept per output coordinate. Its ridge is selected on outer-training people; each test person's two calibration repetitions are explicitly disclosed. This is a deployable calibration scenario, distinct from an unseen-person motion waveform claim.

M4b uses the F0/F1/F2/F3/R descriptor families at the **trial/source-set level** and task indicators to predict `sqrt(meanBandPower)`. Training-only positive medians define asinh feature/target scales, and all fitted transforms exclude test people. A fixed-ridge tonic-energy/MVC comparator is common to the compared methods. Missing method features fall back to that comparator, with missingness retained. Physical squared errors remain in mm² or rad²; transformed losses are dimensionless. A continuous power integration guard of [2.5,9.5] Hz avoids dropping half of a smoothed 3- or 9-Hz line; its exact edges are saved as `powerIntegrationBand`. Search-band detection remains 3–9 Hz.

M5 uses disjoint multitaper windows. Default attribution requires full input rank, condition number <=1e6, nominal replicates >=3 per input, at least five independent windows, and multiple coherence >=0.5. Fifteen correlated EMG inputs in 30 seconds will often fail these conditions. Ordinary coherence remains a descriptive association; regularized multiple coherence is labelled when attribution is unavailable. Contributions include interference, may be negative or exceed one, and are not constrained to sum to one. They describe a fitted linear spectral model, not intervention-based muscle causality. Recorded power quantities here are **one-sided spectral densities at a frequency**, not integrated power or watts.

Surrogates default to zero, so p-values are unavailable. Optional whole-window pairing preserves each taper bundle and local spectra and uses a maximum over searched frequencies/inputs. It requires exchangeability; a stable common line can remain coherent after pairing. It is not a universal test of independent oscillators, causal direction or neural source origin. Sensitivity to window length and input subset belongs in explicitly labelled reruns.

## Selection, severity and operational information loss

All observations from a person remain in one clinical outer fold. Inner held-out people select each EMG configuration and mechanical readout ridge. Existing control/severe/coverage/electrical-calibration gates remain; lower-severity capture comes first. When every clinically tied candidate has finite mechanical validation utility, RMS gain is the optional tie-break before balanced clinical log loss. Turning `cfg.motion.selectionTieBreak=false` performs an EMG-selection sensitivity run. No outer-test mechanical result updates the candidate pool or reference policy.

TETRAS severity definitions and ET/PD phenotype metadata remain the existing explicit clinical fields. Without provided grades or fixed cuts, lower/higher severity means outer-training-relative tertiles. Mechanical tables retain those definitions, reference evaluability and control/ET/lower/moderate/higher strata. The motion-only clinical comparator has its own training-controlled threshold and the same coverage/severe safeguards. It combines only one declared coordinate kind; it does not mix physical units or add motion to an EMG classifier.

For compression, compare the rich descriptor readout R with retrained F0/F1/F2/F3 and with a frozen R readout whose omitted features are replaced by training medians. For paired held-out targets, let B be baseline loss, R rich loss, and C compressed loss. Operational retention is `(B-C)/(B-R)` only when the rich gain exceeds a declared floor. Loss increase is `C-R`. Negative retention and values above one remain visible. Weak/negative rich gain makes the ratio unavailable. Subject bootstrap intervals are conditional on the already fitted outer models; they do not repeat the full nested optimization.

`motion_waveformPrediction.csv` also compares native waveform utility with trial-constant F0/F1/F2/F3/R descriptor readouts and a frozen native readout fed the waveform's temporal mean. `RsubjectWaveform*` denotes the rich **descriptor-family** comparator applied per trial, not a waveform representation. These declared compressions are not inverses of arbitrary statistics. Trial-constant features generally cannot preserve burst timing or phase even if global variance and PSD are preserved.

These percentages quantify utility for a specified target and readout capacity, not total Shannon information in an sEMG recording. Native EMG power remains an electrical estimate; mechanical amplitude does not identify true single-muscle electrical tremor power. Keep raw data, selected waveforms, masks, frequency tracks, events and evidence for future endpoints.

## Storage, reproducibility and outputs

Each run retains its automatic run number/ID, candidate/source hashes, complete parameters and selected method labels. References are frozen once per subject/task/trial under `motion/references`. Compact agreement/features/readout errors are cached per candidate. Only one task's three repetitions of temporary method waveforms are retained during extraction, bounded by candidate batch size. Selected full outputs remain in the existing `selected` directory; all-candidate waveforms are not duplicated on disk.

Motion reference/source/clock/policy changes invalidate mechanical cache acceptance. EMG algorithm hooks receive stripped configurations and keep their original input contracts. Missing motion or uncertain states never remove an otherwise eligible person from EMG clinical evaluation. A required run aborts if no matching canonical motion trial imports at all.

| Output | Contents |
|---|---|
| motion_referenceAudit.csv, motion_coverage.csv | Every planned reference trial/target, missingness, units, states, clock and limit calibration |
| motion_agreement.csv, motion_agreementSummary.csv | M1–M3/M6 by reference witness; subject-first task/severity summaries with evaluable counts |
| motion_waveformPrediction.csv, motion_waveformSummary.csv | Person-calibrated M4a and waveform-compression scope/coverage |
| motion_predictions.csv, motion_summary.csv | Outer subject-held-out M4b physical values, losses and gains |
| motion_retention.csv | Frozen/retrained descriptor retention, severity strata, conditional intervals |
| motion_contrasts.csv | Paired subject differences in F1 mechanical gain |
| motion_residualRepeatability.csv | Within-task repetition SD of transformed prediction error conditional on measured motion |
| motion_policies.csv, motion_waveformPolicies.csv | Training membership, chosen ridges and distinct readout scopes |
| motion_clinicalPredictions.csv, motion_clinicalScorecard.csv | Motion-only clinical comparator |
| motion_coherence.csv | M5 descriptive and conditional-identifiability detail, when enabled |
| motion_calibration.csv | Known-input synthetic motion challenges |
| motion/models, motion/clinical_models | Fitted outer RMS and motion-only clinical models |
| motion/references | Signed physical coordinates, filtered waveforms, states, intervals and independent witnesses |

Trial rows permit paired reanalysis. Cohort summaries average repetitions within task, balance P/K, then average subjects. M4 readout utility also balances named targets without summing their physical units. Unsupported endpoints, missing reference categories and absent strata remain explicit; do not replace NA with zero. Full source-set anatomical identifiability and hardware calibration cannot be inferred from the provided EMG file hierarchy.

## Verification status

Run `START_HERE('validate')` in MATLAB before patient interpretation. Native motion tests include clock drift/phase gating, coordinate units/gaps, pose validity, policy identity, reference calibration gates, guarded resampling, abstention denominators, rank rejection, test-target invariance, the ninth-argument selection regression, waveform readouts, partial merges, caching and an optional synthetic end-to-end study.

The authoring environment has no MATLAB/Octave runtime. Executed checks are independent Python equation/source audits, not execution of MATLAB. `VALIDATION.md` and `validation/motion_numerical_checks.json` state the exact scope. No actual motion, EMG, MVC or TETRAS observations from your cohort were supplied for this update.
