# Module contracts and method availability

## Registered pipelines

`tbench.pipelineRegistry('recommended')` returns the six primary matched pipelines.
`tbench.pipelineRegistry('implemented')` returns all 58 implemented research entries.
`tbench.pipelineRegistry('all')` includes unavailable entries for inspection; preflight refuses to benchmark them.

Each entry contains a stable pipeline ID, analysis ID/family, preprocessing ID, input kind, function handles, explicit parameter structures, expected dynamic capacity, implementation hash and optional external dependency files.

Preprocessing representations are envelope time, envelope spectrum, local spectrum, complex time-frequency, carrier time and selected explicit adapter types. `tbench.study.validatePair` rejects a representation that the named analyzer cannot consume. Default representation/analysis matching is not a proof of globally optimal preprocessing; it is a theory-constrained starting set to compare inside the search.

Several analyses share the same optimal front-end. Detector-domain filtering, identical transforms and fixed reporting observers are memoized. Method labels and parameter hashes remain separate.

## One-muscle boundary

The preprocessor/analyzer receives one channel with its clock, original raw samples, the optional reduced raw branch, reference envelope/energy, validity masks and its own MVC. It does not receive other-muscle waveforms, group labels or TETRAS. Other muscles enter B5 corroboration/readout only after per-muscle inference. Software 4 also removes motion/runtime and clinical/tuning configuration fields from all study EMG hooks through `tbench.motion.inferenceConfig`; custom hooks must use their explicit pipeline parameters. Motion references and joins are separate benchmark modules described in `MOTION_CAPTURE.md`.

A per-record estimator can estimate its own background or frequency. Cohort-trained quantities belong inside the grouped training folds. No source grouping/readout may use held-out clinical answers.

## Output

Start from `tbench.emptyOutput(c)` and satisfy `tbench.validateOutput`:

| Field / capability | Requirement |
|---|---|
| time, valid | Complete shared grid and explicit mask |
| detected | Accepted rhythmic-component flag, not an ET diagnosis |
| frequency, powerMean | Positive/NaN frequency; nonnegative/NaN reported power |
| signedExcessPower | Preserve signed estimates when supplied |
| dynamic | Only true if native proposal supports episodes |
| waveform | Only true if a reconstruction is provided |
| active, frequencyTrack | Required for dynamic output; never active on invalid samples |
| signal, fundamental | Finite on valid support; zero in valid off periods |
| onProbability | Optional model-specific score, not automatically calibrated confidence |
| status, message | Explicit error/insufficient/ok state |
| audit | Native evidence, estimator definitions and provenance |

The study analysis wrapper applies the declared common reporting adapter while retaining internal estimates. Global static methods do not acquire native onset/offset merely because a separate adapter exists.

## Add or compare modules

Use `examples/paired_preprocessing.m` to compare two frozen preprocessing choices with the same analyzer. Give them distinct safe pipeline IDs. `tbench.crossPipelines` and `examples/cross_modules.m` create a compatible preprocessing-by-analysis matrix with frozen transform ownership. Transform kind, single/multiresolution container and Thomson DPSS requirements are checked. Rejected pairs retain a reason. Use `examples/custom_search_space.m` for matched candidate spaces.

A custom hook must be a named function in a source file. Captured anonymous/nested state is rejected: hidden parameters would defeat parameter labels, reproducibility and cache invalidation. Place explicit choices in the parameter structures. Declare source files of external helper dependencies in `pipeline.dependencyFiles`.

Update expected dynamic capacity when replacing an analyzer. Preserve the output contract and correct power reporting adapter. If adding a wholly new catalog analysis, also implement its dispatcher path, input contract and capacity declaration; a wrapper name alone is not an implementation. Restart/clear the catalog after changing its manifest.

Signal and benchmark hook source hashes, external dependency hashes, the package source hash, raw/MVC content and processing parameters enter cache identity. Readout ridge is tuned separately. Clinical labels/readouts are never reused from signal caches; a label change triggers fresh fold/readout fitting.

## Catalog correction

The broader saved prototype had 112 wrapper names, of which multiple family backends were missing. The corrected manifest reports:

| Readiness | Count | Automated study behavior |
|---|---:|---|
| implemented | 58 | Available for native validation and prespecified evaluation |
| not_implemented | 46 | Rejected before expensive work |
| requires_model | 5 | Rejected; trained weights and the missing learned backend are required |
| unsupported_input | 3 | Rejected for the present single-channel input |

See `method_coverage.csv` for every ID, description, review section and status. Unimplemented packet/latent/oscillator/decomposition/combination wrappers are not scored as working SOTA algorithms. Source implementations in the available catalog are compact research variants; they require comparison with reference implementations before attributing a result to an entire published family.

The default six include the existing conventional comparator and interpretable dynamic alternatives. Declaring all 58 primary would increase selection/multiplicity and compute costs; retain the short primary comparison and treat extensions as exploratory.

## Memory, caching and storage

The extraction loop loads one participant at a time within batches of at most 32 candidates and writes per-subject/candidate feature caches atomically. Batch size bounds retained feature structures and may cause extra file reads. It holds small subject-level views for clinical tuning rather than every candidate waveform for the entire cohort. Refinement extracts training participants first; selected held-out features are filled after selection. A complete candidate set still costs computation and disk; reduction does not make the study free.

MUSIC, ESPRIT and matrix pencil pool delay covariance over all supported contiguous runs, retaining original timestamps and never joining across gaps. Native inference coverage is audited. Atomic norm is an explicitly labelled overlapping-block finite-iteration proposal covering the supported record; it is not a full-record exact SDP. Static outputs still have no native episodes. IRASA's intentional spectral warp uses anti-alias filtering/guarded sinc sampling and remains separate from the physical reporting-domain reducer.

Only selected outer-test waveforms are retained by default. Native time axes are stored once, fixed observer arrays once per trial, and method outputs by run/fold/pipeline/candidate/subject. Load these through `tbench.loadSelected`. Original inputs and per-record estimates remain available for additional investigations.

Cache hits require content/config/source identity, not filenames or modification times alone. Corrupt/incomplete cache payloads are recomputed. A fresh run recomputes clinical fits and produces a new run number. Fold checkpoint files record progress; the runner does not pretend to resume a partly fitted clinical run in place.

Timings inside extraction reflect transform sharing and cache context. They are workflow timing diagnostics, not controlled independent speed comparisons of algorithms.

Disable optional explanatory work or selected waveform storage through `cfg.execution` when appropriate. Turning those off limits subsequent signal/physiological audits; it does not change the primary signal extraction, grouped clinical comparison or calibration definitions.

