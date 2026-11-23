% A prespecified search replacement with parameter coverage and refinement.
cfg=tbench.studyConfig(); cfg.tuning.maxCandidates=32;
space=struct('analysisId','hann_tracker','preprocess',struct('rawHighpassHz',[20 35], ...
    'localWindowSeconds',[1.5 2 3]),'analysis',struct('threshold',[2 3 4 6], ...
    'minCycles',[3 4],'switchPenalty',[.5 1 2],'frequencyPenalty',[.15 .35 .7]));
cfg.tuning.spaces=space; cfg.tuning.ridgeValues=[1 10 100];
% All choices must be declared before looking at outer-test results.
% For a small complete Cartesian search instead, set searchMode='grid' and
% raise maxGridCandidates explicitly if needed. Coverage is the bounded default.
