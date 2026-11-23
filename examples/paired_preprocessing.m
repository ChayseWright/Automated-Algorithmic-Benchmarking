% Declare a second PREPROCESSING choice without changing the analysis method.
cfg=tbench.studyConfig(); p=tbench.pipelineRegistry('hann_tracker');
p2=p; p2.preprocessId='local_hann_highpass35'; p2.preprocessParams.rawHighpassHz=35;
p2.id=[p2.preprocessId '__' p2.analysisId];
cfg.pipelines=[p p2]; cfg.tuning.enabled=false; % Fixed comparison, no cohort selection.
% The common reporting envelope stays at the study's fixed 20--450-Hz carrier.
% study=run_benchmark('/absolute/data','/absolute/clinical.csv',cfg);
