% RUN_COMPARISON Edit the folder and MVC manifest, then execute this script.
% EMGData.P4.trial2 and EMGdata.P4.trial2 are both supported.
addpath(fileparts(fileparts(mfilename('fullpath'))));
cfg=tbench.defaultConfig();
subjects=tbench.subjectsFromFolder('CHANGE_ME_TO_YOUR_SUBJECT_MAT_FOLDER');
% If each file contains MVC (15-element EMG amplitude vector), it is loaded.
% Otherwise assign references, in the SAME column order as the 15 muscles:
% subjects(1).mvc = MVC_by_subject(1,:);
% cfg.muscleNames = {'your first muscle', ... };  % Exactly 15 names.
% cfg.nominalFs = 1111;                          % Optional timestamp check.
% cfg.timeScale = 0.001;                        % ONLY if time is milliseconds.
cfg.synchronized=true;  % Verify that the 15 columns share the same acquisition clock.
% preview=tbench.inventory(subjects,cfg);
study=tbench.run(subjects,cfg);
