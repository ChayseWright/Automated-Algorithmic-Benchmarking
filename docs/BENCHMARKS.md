# Implemented benchmark definitions and capacity

The unchanged source is `tremor_semg_benchmark_protocol.pdf`, revision November 3, 2026, SHA-256:

`0fa4b3fb0abaf54645988cb65a991d18cda849ba4c9ea05e9b1debe1d35c5be7`

The current user request overrides two earlier defaults: the primary search is **3–9 Hz**, and **lower-severity capture** is the first selection objective. The 4–8 Hz endpoint is secondary; 2–15 Hz is a fixed diagnostic witness. The benchmark is an empirical research comparison, not ET-versus-PD diagnosis from a frequency cutoff.

## Acquisition hierarchy and aggregation

The inference unit is one participant. Subject-level features pool eligible trial–muscle cells within each task, average tasks within postural/kinetic types, and give those two types equal weight. Missing task/record counts and actual valid fractions remain visible. A participant must meet the declared coverage requirement and have both task types to supply a definitive primary readout.

A record's accepted-rhythm flag is O. The clinical subject decision D comes from a separate small readout and training-fitted threshold. An ET participant is not automatically positive in all 585 muscle records, and a control participant is not automatically free of physiological rhythm.

## G0 — Measurement eligibility

Connected modules: `validateTrial`, `preprocess.reference`, `study.metrology`, `study.reductionAudit`.

Checks include timestamp/sample-rate consistency, missing-sample expansion, declared-rail masks, flat channels, finite MVC references, output capabilities, gain/polarity behavior, actual anti-alias response and support. The fixed reference is computed before carrier reduction. Actual sampling rates, output rates, valid fractions and runtime are retained.

Changing raw gain by c should multiply electrical power by c²; scaling MVC by |c| should preserve MVC-normalized power. Metrology retains the gain ratio, frequency/detection changes and response diagnostics. The selected-output reference preserves the original time origin; episode times are on the relative trial clock.

No automatic mask is inferred from “low TETRAS,” high amplitude, or apparent diagnostic improvement. Undeclared clipping, poor electrode placement and unreliable MVC cannot be perfectly established from the given data alone. Review masks are retained when an algorithm explicitly supplies them; a zero review fraction is not a validated clean-signal label.

## B1 — Controls

Primary endpoint: number of **definitive subject-level ET-like positives / all enrolled controls**. Also report severe-label calls in controls, control coverage and indeterminate counts.

Secondary control physiology: task-specific and task-balanced accepted rhythm yield in the primary band and 4–8 Hz, valid-time/MVC power, occupancy and episode rate where supported. `rhythmSummary.csv` includes individual participants as the aggregation units; do not treat 5,850 control muscle records as independent subjects.

Control rhythms remain measurements, not automatically DSP errors. Generated no-modulation alarms in B4 have a different reference meaning from natural control rhythms.

The threshold is selected only from inner out-of-fold training-control scores. The requested 0.10 operating point is empirical. With 10 actual controls, outer FPR changes in 0.10 steps; repeated folds cannot establish rare-error specificity. No exact-binomial population guarantee is attached to the dependent cross-validated predictions.

## B2 — Lower/higher-severity capture

“Mild capture” = definitive ET-like positives among independently labelled lower-severity ET participants / **all** such participants. An indeterminate output is not a captured case. Higher-severity capture is defined analogously. These are case-capture endpoints, not sensitivity for every natural EMG episode.

`clinicalScorecard.csv` retains numerators, denominators, coverage and indeterminate counts. `operatingCurves.csv` gives descriptive outer-score ROC/capture points. `severityCaptureCurves.csv` reports cumulative capture below/above actual severity-score values at the selected outer operating points. `withinFoldAUC.csv` checks sensitivity to pooling scores from different fits.

Provided clinical categories or fixed, justified cuts are preferred. Otherwise relative tertile cutpoints are fitted separately to training ET participants. Ties remain together; strata can be unequal or empty. The exact score kind and fold model/cutpoint provenance remain available. Do not reinterpret a relative low tertile as a universal clinical definition of mild ET.

The default three-variable common readout uses task-balanced mean power, upper-record power and accepted-record fraction. Upper-record power is a within-task maximum by default, retaining a rare affected channel/trial. It can be influenced by artifact; the control/fidelity panel must accompany it.

## B3 — Clinical severity and asymmetric error

Within ET only, the suite reports Spearman correlation, Kendall tau-b, numeric-score MAE where enough scores exist, conditional ordinal ranked probability score, and severe undercalling/false-severe counts/rates.

Higher-severity “capture” is any positive ET-like case flag. It is **not** the same as correctly assigning the severe grade. A severe case labelled mild is captured but undercalled. Indeterminate severe cases are reported separately.

The four-class cost matrix indexes true/decision classes [control, lower, middle, higher]. The default penalizes missing a higher-severity case more heavily than missing a lower-severity case. The case flag remains constrained by the separately chosen control operating point; severity labels for positive flags minimize the declared conditional decision cost.

The primary binary selection objective does not optimize every severity endpoint simultaneously. False-severe calls and within-ET agreement remain explicit axes. The costs/cutoffs are prespecified research choices, not numerical clinical facts learned from this cohort.

The tonic/MVC baseline (B0 below) helps identify nuisance separation. Age, device, medication, side and task matching require prespecified metadata policies and additional analyses; they are preserved, not automatically adjusted away. A global subject score cannot substitute for a missing same-arm/task TETRAS annotation.

## B4 — Known-input accuracy

`study.calibrate` operates on generated raw Gaussian-carrier AM signals, not directly injected envelope sinusoids. Its target is the **artifact-free expected reporting-envelope projection**, derived through the fixed original-rate carrier filter, rectification expectation and anti-alias observer.

The default 75 cases cross three frequencies, five depths and five temporal patterns:
- 3, 6 and 9 Hz;
- no modulation, 0.05, 0.10, 0.25 and 0.80 depth;
- continuous, early, middle, late and repeated short episodes, with the shortest gates retaining at least three cycles and at least 0.5 s;
- selected frequency drift, transient contamination, missing intervals and declared clipping.

Synthetic amplitude units are independent of the participant ADC scale. Participant rail thresholds are disabled for generated data; designated clipping challenges use their own declared rails. A clipped/gapped challenge may legitimately abstain.

Development cases gate candidate selection; separate seeds challenge selected outer configurations. Rows retain detection, status, known/estimated mean power, relative error, frequency error, occupancy confusion counts, episode IoU, onset/offset error, conditional track error and power on the truth-active samples. Missing frequency estimates are not counted as exact estimates. Static methods return NA for absent native timing.

“Conditional track RMSE” is conditional on a detected true-active overlap and must be read beside missed/false-active samples. Onset/offset errors summarize the first/last accepted support. For repeated bursts these are record-span boundary errors, not matched-event boundary accuracy; false-positive samples and IoU remain necessary. Static methods have no native episode score. Separate weak-depth frequency-recovery and late strong-burst recovery fractions expose failures that an overall average can hide.

The generator is a limited development challenge set, not an exhaustive model of motion artifact, every motor-unit physiology or every ET/PD phenotype.

### Optional held-out control counterfactual

Enable `cfg.bench.controlInjection=true`. Default probes use P1/K1, trial1, technical columns 1/8/15, frequencies 3/6/9 Hz, and depth 0.10/0.80. These channel indices have no implied anatomical selection.

For an outer-test control, multiply its original raw carrier by a known gated modulation, then recompute the complete observer and method. Truth is the known-phase fundamental projection of the **exact counterfactual reference-envelope change**. This deliberately includes nonlinear interactions with that control background and is not identified latent neural drive.

Where both method outputs contain waveforms, project their waveform difference onto the same known phase to compare incremental power. Static methods cannot supply that comparison; their incremental-waveform result is NA. A difference in their total power is retained descriptively and is not passed off as incremental component power.

Control backgrounds are not reused to tune the configuration being tested. Several probes in one control remain dependent measurements from one participant.

## B5 — Physiological support

Connected checks include:
- Held-out target-energy predictive gain from other-muscle outputs, with training-only source selection/readout.
- Reference/candidate complex coherence and phase relationships.
- Independent fixed envelope/energy frequency witnesses from other muscles.
- Residual-band power and the exact variance cross-term identity.

Cross-muscle prediction/coherence require `cfg.synchronized=true`; do not assert it solely because columns share a matrix. Frequency witnesses and residual checks can still be reported. Corroboration is not a physiological ground truth: common artifact, cross-talk and task drive can generate support in both ET and controls.

Blocked readout training is guarded by at least the actual fixed energy-observer filter support. Source inference uses its own muscle; the held-out target does not enter source estimation. This is concurrent cross-muscle prediction, not a claim of forecasting an unseen future episode.

Optional `cfg.bench.temporal=true` uses raw-target masking before every filtering/inference step, followed by a named **masked-frequency static readout**. It tests retained frequency information, not native burst forecasting. Baseline fallback and unavailable prediction are retained rather than scored as a perfect reconstruction.

Detailed B5 rows stream to CSV. Subject-first, task-type-balanced summaries and participant paired comparisons are separate from pooled within-person samples.

## B6 — Repeatability

Within participant/task/muscle, report repetition detection discordance, paired frequency differences and positive-power log differences. Zero and missing powers are not converted to arbitrary log-epsilon values.

ICC(A,1), absolute SEM and MDC95 are calculated for complete repeated data separately by task/muscle and diagnostic group, for absolute and MVC power. MDC95 = 1.96*sqrt(2)*SEM. The SEM uses residual variance plus a nonnegative occasion-variance estimate. At least three complete participants and two repetitions are required.

MVC-based total power and 15-muscle composition changes are retained separately. Composition requires complete nonnegative power and positive totals; a similar topography does not imply similar total burden.

Outer participants can have different training-selected parameter sets. These B6 results describe the selected workflow across participants, not an ICC of one universally frozen parameter vector. A fixed-parameter sensitivity run is available through disabling tuning.

## B7 — Processing and feature compression

Two related but distinct ladders are implemented.

### Clinical feature ladder

| Level | Features | Width / capacity |
|---|---|---|
| B0 | Tonic energy and MVC amplitude | 2; signal-independent nuisance comparator |
| F0 | Envelope variance and accepted-record fraction | 2; original benchmark comparison |
| F1 | Mean rhythmic power, upper-record power, accepted-record fraction | 3; common primary |
| F2 | F1 + occupancy, upper active-window power, activity indicator | 6; dynamic only |
| F3 | F2 + episode power, duration, episode rate | 9; dynamic only |
| R | F3 + local-power median/upper tail, harmonic ratio, frequency spread/change, envelope variance | 15; richest supplied subject summary |

Within a record, active-window power uses one-second windows, or the whole accepted episode when it is shorter. Episode power is mean fundamental-square within a contiguous accepted episode. Episode energy includes power times duration. Per-record maxima, quantiles and all retained episode rows preserve additional detail outside the small clinical feature vector.

The zero used for an off-state descriptor is a **model encoding**, accompanied by an activity indicator. It is not an assertion that an undefined active power is a measured zero. A static R contains its supported static features/variance and constant placeholders for absent descriptors. It cannot recover missing episode information.

F0 defines variance in the fixed reference envelope and retains variance/yield as separate readout inputs. It does not assume an unspecified product formula, or fabricate a reconstructed-waveform variance for a static method.

Two clinical ablations are exported:
1. **Training-recalibrated:** refit the same class of clinical readout to each feature level using training participants only, and choose its threshold from training controls.
2. **fixed_F0–fixed_F3:** keep the rich model and threshold fixed; replace omitted raw descriptors by training medians, retaining matched descriptor values and masks.

Signal parameters/ridge are selected on the primary level, not separately optimized using each outer feature-ablation result. This is a conditional audit of the selected signal workflow. B0 is an alternative nuisance baseline rather than a nested compression of R.

For binary balanced log loss L, compare a baseline ln(2), rich loss LR and compressed loss LC. Loss increase is LC−LR. Retained diagnostic utility is (ln(2)−LC)/(ln(2)−LR) only if rich gain is positive. Do not clip it: negative values or values above one are possible. Mild capture, control FPR and coverage changes accompany it.

### Physiological ladder

B5 compares yield-only, static frequency/power, episode summaries and rich waveform/quadrature predictors at equal predictor width where supported. It reports held-out predictive-gain loss and retained utility relative to its defined rich representation.

Clinical R and physiological “rich” are different references: one is a subject descriptor vector; the other retains a waveform. Neither is all information in raw EMG.

These are **operational, target/readout-specific information audits**. Finite-sample decoder limitations, regularization, missingness and target choice affect them. A single scalar variance/yield benchmark cannot determine a universal percentage of lost Shannon information, and two signals with identical benchmark features can have different phase/timing/coherence. Original recordings and selected continuous outputs remain the route to testing additional targets.

## B8 — Applicability and phenotype

The suite retains method errors, insufficient support, indeterminate clinical decisions, group-specific coverage, task-level failure counts and fixed broad-band witnesses outside the 3–9 Hz primary search. It reports all-enrolled results and a secondary common-definitive-output participant comparison.

Phenotype scorecards evaluate each documented ET phenotype plus controls under the **inclusive training policy**. They do not silently retrain after exclusions. To evaluate a restricted ET/ET-plus policy, freeze that manifest and run a separate grouped experiment. PD-only or uncertain-diagnosis populations require an explicitly extended label/evaluation design; do not relabel them from a PSD peak.

Neither ET–PD attribution nor ET-versus-PD specificity is established by this dataset without the corresponding independent clinical references. No prevalence of ET–PD co-occurrence is assumed.
