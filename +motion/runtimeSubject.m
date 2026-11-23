function [files,key] = runtimeSubject(cfg,id)
%RUNTIMESUBJECT ID-based lookup also works for outer-training subsets.
files=cell(numel(cfg.tasks),numel(cfg.trials)); key='motion_disabled';
if ~cfg.motion.enabled, return; end
if ~isfield(cfg,'motionRuntime')
    error('tbench:MotionRuntime','Build independent motion references through runStudy before motion-enabled extraction.');
end
r=cfg.motionRuntime; s=find(strcmp(r.subjectIds,id),1);
if isempty(s), error('tbench:MotionJoin','Subject missing from frozen motion manifest.'); end
for j=1:numel(cfg.tasks), for k=1:numel(cfg.trials), files{j,k}=r.files{s,j,k}; end; end
key=tbench.study.hash(struct('motionSource',r.subjectHashes{s},'policy',cfg.motion));
end
