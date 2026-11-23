# Parameter selection and experiment discipline

## What “best” means

A configuration consists of a named preprocessing recipe, analysis ID, signal parameters, clinical feature level and ridge penalty. Candidate hashes identify signal settings and module implementations. Configuration hashes additionally identify clinical settings and ridge values. Run/fold identifiers identify the actual fitted experiment.

There is no data-independent universal optimum. The software finds the best **evaluated** configuration under the declared benchmark. Software 3 exhausts small discrete spaces; larger spaces begin with the prespecified baseline and one-coordinate anchors covering every declared value, then reproducibly fill a 64-configuration budget. Two training-only refinement rounds propose up to 16 configurations each around up to three leading configurations. These rounds explore interactions without guaranteeing a global optimum or convergence.

Do not tune the benchmark to make a method win. Freeze the clinical target, eligibility masks, severity definition, control budget, aggregation, candidate families, search spaces and cost policy before examining outer results.

## Default spaces

| Parameter | Applicable branch | Candidate choices |
|---|---|---|
| Raw high-pass | Detector domain | 20, 35, 40 Hz |
| Welch window | Static spectrum | 2, 3, 4, 6 s |
| Local window | Local/seeded tracking | 1, 1.5, 2, 3, 4 s, where applicable |
| Local hop | Dynamic proposals | 0.1, 0.2, 0.3 s |
| Evidence ratio | Applicable detectors | 2, 3, 4, 6, 8 |
| Minimum cycles | Dynamic proposals | 3, 4, 5 |
| Switching penalty | Local tracking | 0.5, 1, 2 |
| Frequency penalty | Local tracking | 0.15, 0.35, 0.7 |
| Wavelet cycles | Applicable CWT/BOSC variants | 3, 4, 5, 7 |
| Morse gamma | Morse/BOSC | 2, 3, 4 |
| DPSS bandwidth / tapers | DPSS, multitaper, Thomson, Kalman seed | NW 1.5/2/2.5/3; K 1/2/3/4, repaired to K <= 2NW-1 |
| Spectrum overlap | Applicable spectral methods | 0, 0.5, 0.75 |
| Multiresolution windows | Multiresolution STFT | Long 2/3/4 s; short 0.75/1/1.5 s |
| Kalman stay probabilities | Switching model | 0.98, 0.995, 0.999 |
| Process fraction | Switching model | 0.0003, 0.001, 0.003 |
| Student innovation parameter | Switching model | 3, 5, 10 |
| Observation fraction / damping | Switching model | 0.05/0.15/0.30; 0.995/0.999/0.9995 |
| AR order | AR methods | 4, 8, 12, 20 |
| Subspace embedding / rank | MUSIC, ESPRIT, pencil | 24/40/64/96; 2/4/6 |
| Atomic block / penalty / iterations | Block atomic proposal | 128/192/256; 0.1/0.3/0.5/0.8; 40/80/160 |
| Guard / reference / background quantile | Applicable background estimators | 0.5/0.75/1 Hz; 3/5/7 Hz; 0.4/0.5/0.7 |
| Clinical ridge | Every signal candidate | 1, 10, 100 |

The authoritative method-by-method definition is `+tbench/+study/parameterSpace.m`, also exported as `parameter_spaces.csv`. Spectrum, AR, subspace, transform, cyclic and demodulation families expose different parameters. Unknown custom analysis fields are rejected. Variant-specific spaces may be replaced through `cfg.tuning.spaces`; inspect the implementation before extending the supported field list. The catalog does not claim optimization over all possible algorithms, priors or continuous parameter values.

The reporting envelope, analysis clock, harmonic count and physical power definition remain frozen. Candidate bandwidth/support are checked. DPSS pairs are repaired to their smallest admissible NW; applicable local windows are enlarged to fit the declared minimum cycles at 3 Hz. Morse beta = cycles^2/gamma. Repaired candidates are deduplicated by their final content hash. Cross-module experiments freeze the chosen preprocessor's transform parameters rather than searching transform parameters that the analyzer would ignore.

`cfg.tuning.searchMode='grid'` requests full Cartesian enumeration, guarded by `maxGridCandidates=4096`. `requireParameterCoverage=true` rejects an insufficient coverage budget. Reduced examples explicitly turn that requirement off and disable refinement. An early stop means no unseen coordinate neighbor within this search scope; hitting the refinement budget is reported as such, not described as convergence.

Use `examples/custom_search_space.m` to replace a method's space. Signal parameters must be explicit, not captured in anonymous functions. Custom external dependencies must be listed for hashing.

## Nested execution

1. Validate metadata, module availability and representation contracts.
2. Run native installation probes and development known-input calibration.
3. Extract all prespecified label-free signal candidates once. Cache them by source content, MVC, code/module signatures, processing configuration and parameters.
4. Make outer participant folds, stratified by group and distributed over ordered severity.
5. For each outer training set, make inner participant folds. Learn power transforms, imputation/standardization, clinical models and severity cuts on inner training participants only.
6. Generate inner out-of-fold scores and choose a strict `score > threshold` rule from eligible training controls.
7. Rank initial configurations using inner results and development calibration, then refine only within that outer-training set. Refinement candidates proposed by another fold never enlarge this fold's eligible search pool. Identical label-free physical calculations may reuse caches.
8. Extract a newly selected candidate on outer-test participants only after selection; refit on outer training participants and evaluate outer participants once.
9. Run independently seeded calibration on the chosen outer configurations.
10. Fit a separate all-training policy for later use.

With 25 ET/10 controls and five folds, each test fold contains five ET and two controls. Each training fold has 20 ET/eight controls. Four inner folds are used when counts permit. Participant IDs and every inner train/validation membership are exported.

The clinical model is a small ridge logistic readout with a balanced working prior. Its scores are not automatically calibrated probabilities at clinical prevalence. Severity uses conditional ordinal models and an optional numeric score readout. The numeric readout stays within its training score range.

## Default acceptance gates

These are engineering development targets, not established clinical standards:

- Empirical inner control-positive rate no greater than 0.10.
- At least 0.80 clinical coverage overall and in each diagnostic group.
- At least 0.80 higher-severity case capture.
- A lower-severity stratum must be evaluable when it is the selection priority.
- Known-input status coverage at least 0.90.
- Strong intermittent challenge detection at least 0.70.
- Median strong intermittent relative mean-power error at most 0.75.
- Median strong intermittent frequency error at most 0.75 Hz.
- Record-level accepted flags in generated no-modulation challenges at most 0.50.

The 75-case calibration grid is a finite development set. It crosses three frequencies, five depths and five episode patterns; designated cases also include artifacts. Show case counts and failures; these thresholds do not imply a confidence bound. Weak-depth and late-burst recovery are reported separately. Development calibration informs selection; independent-seed selected-configuration calibration is never fed back into selection.

Among eligible configurations, maximize lower-severity capture. When motion is enabled and `selectionTieBreak=true`, tied candidates use inner held-out RMS-prediction gain only if every tied candidate has a finite comparable motion result. Remaining ties use lower balanced log loss and then original candidate order. `cfg.clinical.priority='severe'` restores the PDF's earlier high-severity-first objective. The frozen reference policy itself is not tuned on test people; see `MOTION_CAPTURE.md` for distinct population and person-calibrated readout scopes.

If no candidate is eligible, minimize the recorded aggregate constraint violation, then apply the same tie policy, and label the result `constraints_unmet`. This fallback is an auditable provisional configuration, not a passing winner or a silent relaxation. Calibration errors and missing predictions remain in the tables.

The default upper-record burden is a maximum within each task before task-type balancing (`upperRecordQuantile=1`). This keeps a rare affected muscle/trial visible. It also exposes artifact/noise extremes, so control outcomes, known challenges and quality checks remain necessary. A 0.90 upper quantile is a prespecified sensitivity option; choosing it after seeing outer outcomes would invalidate the primary selection claim.

## Control-threshold resolution

For n eligible inner out-of-fold controls, allow floor(alpha*n) empirical positives; choose the corresponding descending control score and use strict exceedance. Ties cannot create extra exceedances.

At alpha = 0.10 and eight controls, the permitted empirical count is zero. The fitted model can shift after refitting; the outer control FPR is what actually matters. This is not an independent split-conformal threshold, and no population-specificity guarantee is attached. Repeated folds do not add independent controls.

## Reporting and uncertainty

Read `selected_parameters.csv` by its scope column: outer-training selections versus all-training deployment choices. A final all-data parameter vector has not itself been independently tested. The outer results estimate the tuning-and-fitting procedure.

Paired contrasts resample participants within control and ET severity strata, carrying each person's predictions together. Intervals are conditional on the already fitted outer models; they omit full model-selection uncertainty. Within-fold AUC is provided alongside pooled score AUC. Operating/capture curves are descriptive outer-result displays and are not fed back into tuning.

`optimizationTrace.csv` records the training IDs, parent/winner hashes, candidate counts, search round and termination reason. `parameter_stability.csv` describes outer-fold choices; the overlapping training sets are dependent. These stability summaries are not independent replications or a reason to retune the outer benchmark.

Do not select the method with the highest outer result and call its selected maximum an unbiased final performance estimate. Confirm a provisional method choice in an independent cohort/session.

