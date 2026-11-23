function study = run_motion_benchmark(dataFolder,clinicalCSV,motionFolder,cfg)
%RUN_MOTION_BENCHMARK EMG study plus independent matched motion references.
% Empty motionFolder uses embedded MotionData in each subject EMG MAT file.
if nargin<3, motionFolder=''; end
if nargin<4, cfg=tbench.studyConfig(); end
cfg.motion.enabled=true; cfg.motion.required=true;
subjects=tbench.attachClinical(tbench.subjectsFromFolder(dataFolder),clinicalCSV);
subjects=tbench.attachMotion(subjects,motionFolder);
cfg.clinical.sourceFile=clinicalCSV; cfg.clinical.sourceSHA256=tbench.study.fileHash(clinicalCSV);
study=tbench.runStudy(subjects,cfg);
end
