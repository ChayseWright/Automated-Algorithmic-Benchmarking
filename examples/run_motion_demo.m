% Small, explicitly synthetic integration run for the motion extension.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
cfg=tbench.studyConfig(); cfg.expectedSubjects=12; cfg.expectedDuration=16;
cfg.tasks={'P1','K1'}; cfg.trials={'trial1','trial2','trial3'}; cfg.outerFolds=3;
cfg.tuning.innerFolds=2; cfg.tuning.maxCandidates=2; cfg.tuning.requireParameterCoverage=false;
cfg.tuning.refinementRounds=0; cfg.tuning.ridgeValues=10; cfg.tuning.bootstrapSamples=50;
cfg.bench.calibration=false; cfg.tuning.calibrationGate=false; cfg.bench.metrology=false;
cfg.execution.computeExplanatory=false; cfg.motion.enabled=true; cfg.motion.required=true;
cfg.synchronized=true; % This synthetic generator uses one EMG acquisition clock.
cfg.motion.minimumFitRows=8; cfg.motion.minimumMotionSubjects=2;
cfg.motion.readoutRidges=[1 10]; cfg.motion.calibration.enabled=false;
cfg.motion.coherence.minimumSegments=3; cfg.motion.coherence.surrogates=0;
cfg.pipelines=tbench.pipelineRegistry({'welch_cfar','hann_tracker'});
cfg.outputDir=fullfile(tempdir,'tbench_motion_demo'); cfg.cacheDir=fullfile(cfg.outputDir,'cache');
subjects=tbench.demoMotionSubjects(cfg,6,6);
study=tbench.runStudy(subjects,cfg);
disp(study.outputDir);
