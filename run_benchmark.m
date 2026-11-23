function study = run_benchmark(dataFolder,clinicalCSV,cfg)
%RUN_BENCHMARK Main study entry; one MAT and one clinical row per person.
if nargin<2, error('tbench:Usage','Provide the subject-MAT folder and clinical CSV filename.'); end
if nargin<3, cfg=tbench.studyConfig(); end
subjects=tbench.subjectsFromFolder(dataFolder);
subjects=tbench.attachClinical(subjects,clinicalCSV);
cfg.clinical.sourceFile=clinicalCSV; cfg.clinical.sourceSHA256=tbench.study.fileHash(clinicalCSV);
study=tbench.runStudy(subjects,cfg);
end
