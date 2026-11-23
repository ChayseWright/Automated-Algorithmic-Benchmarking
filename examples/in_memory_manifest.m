% IN_MEMORY_MANIFEST Adapt a 25-subject cell array already in your workspace.
% Required inputs: EMG_by_subject{s}.P4.trial2, MVC_by_subject(s,1:15).
projectRoot=fileparts(fileparts(mfilename('fullpath'))); addpath(projectRoot);
subjects=repmat(struct('id','','EMGData',struct(),'mvc',[]),25,1);
for s=1:25
    subjects(s).id=sprintf('S%02d',s);
    subjects(s).EMGData=EMG_by_subject{s};
    subjects(s).mvc=MVC_by_subject(s,:);
end
cfg=tbench.defaultConfig(); cfg.synchronized=true;
study=tbench.run(subjects,cfg);
