function [data,identity] = loadSubject(spec,m)
%LOADSUBJECT Canonical import only; no inferred vendor channel/axis mapping.
data=struct(); identity=struct('status','missing_motion','source','','hash','');
if isfield(spec,'MotionData')&&~isempty(spec.MotionData)
    data=spec.MotionData; identity.source='in_memory_MotionData';
    identity.hash=tbench.study.hash(data);
else
    file=tbench.option(spec,'motionFile','');
    if isempty(file), file=tbench.option(spec,'file',''); end
    if isempty(file)||exist(file,'file')~=2, return; end
    variables=whos('-file',file); names={variables.name};
    variable=tbench.option(spec,'motionVariable',m.variable);
    if ~any(strcmp(names,variable)), return; end
    raw=load(file,variable); data=raw.(variable);
    identity.source=file; identity.hash=tbench.study.fileHash(file);
end
if ~isstruct(data)||~isscalar(data), error('tbench:MotionRoot','MotionData must be one task structure per subject.'); end
identity.status='ok';
end
