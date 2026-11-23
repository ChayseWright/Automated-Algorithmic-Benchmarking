% Compatible preprocessing/analysis experiments, declared before evaluation.
cfg=tbench.studyConfig();
P=tbench.pipelineRegistry({'hann_tracker','dpss_tracker'});
[cfg.pipelines,rejected]=tbench.crossPipelines(P,P,cfg);
% Four pairs: Hann or DPSS evidence consumed by each compatible tracker.
% Transform parameters belong to the selected preprocessor and are frozen.
% Threshold/duration/path settings and carrier high-pass remain searchable.
disp({cfg.pipelines.id}'); disp(rejected);
% study=run_benchmark('/absolute/data','/absolute/clinical.csv',cfg);
