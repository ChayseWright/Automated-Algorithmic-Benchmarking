% RUN_COMPARISON Software-3 clinical/control benchmark, one file per participant.
addpath(fileparts(mfilename('fullpath')));
cfg=tbench.studyConfig();
dataFolder='CHANGE_ME_TO_YOUR_SUBJECT_MAT_FOLDER';
clinicalCSV='CHANGE_ME_TO_YOUR_CLINICAL_CSV';
% Map your exact 15 acquisition columns, time units, and measured MVC references.
% cfg.muscleNames={'first muscle', ...}; % Exactly 15 names.
% cfg.nominalFs=1111;                   % Optional check against the time column.
% cfg.timeScale=.001;                  % ONLY for millisecond timestamps.
% cfg.synchronized=true;               % ONLY after verifying the shared clock.
% cfg.clinical.severityMode='provided'; % When clinician-defined categories exist.
% cfg.clinical.severityCuts=[... ...];  % Alternatively, prespecified score cuts.
subjects=tbench.subjectsFromFolder(dataFolder);
subjects=tbench.attachClinical(subjects,clinicalCSV);
cfg.clinical.sourceFile=clinicalCSV; cfg.clinical.sourceSHA256=tbench.study.fileHash(clinicalCSV);
preview=tbench.inventory(subjects,cfg); % Review missing records/MVCs before a full run.
study=tbench.runStudy(subjects,cfg);
