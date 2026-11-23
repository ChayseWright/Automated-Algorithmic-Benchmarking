% RUN_DEMO A small synthetic smoke run; rows are explicitly labeled DEMO.
projectRoot=fileparts(fileparts(mfilename('fullpath'))); addpath(projectRoot);
cfg=tbench.defaultConfig();
cfg.expectedSubjects=2; cfg.tasks={'P1','K1'}; cfg.synchronized=true;
cfg.bench.calibration=false; % Run the calibration panel separately for a quick demo.
cfg.bootstrapSamples=100;
subjects=tbench.demoSubjects(2,12,1000);
study=tbench.run(subjects,cfg);
disp(study.outputDir);
% To run only the known-generator calibration:
% calibration=tbench.calibration(cfg.methods,cfg);
