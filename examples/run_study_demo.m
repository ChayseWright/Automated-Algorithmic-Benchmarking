% Small synthetic integration run, including nested search and compression.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
cfg=tbench.studyConfig(); cfg.expectedSubjects=12; cfg.expectedDuration=14;
cfg.tasks={'P1','K1'}; cfg.trials={'trial1','trial2','trial3'};
cfg.outerFolds=3; cfg.tuning.innerFolds=2; cfg.tuning.maxCandidates=2;
cfg.tuning.requireParameterCoverage=false; cfg.tuning.refinementRounds=0;
cfg.tuning.ridgeValues=[1 10]; cfg.tuning.bootstrapSamples=100;
cfg.bench.calibration=false; cfg.tuning.calibrationGate=false;
cfg.bench.metrology=false; cfg.execution.computeExplanatory=false;
cfg.clinical.minSubjectCoverage=.50; cfg.outputDir=fullfile(tempdir,'tbench_demo_results');
cfg.cacheDir=fullfile(cfg.outputDir,'cache');
cfg.pipelines=tbench.pipelineRegistry({'welch_cfar','hann_tracker'});
subjects=tbench.demoStudySubjects(cfg,6,6);
study=tbench.runStudy(subjects,cfg);
disp(study.status); disp(study.outputDir);
