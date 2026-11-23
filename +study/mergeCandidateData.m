function data = mergeCandidateData(data,addition,subjectIndices,candidateIndices,cfg)
%MERGECANDIDATEDATA Insert label-free partial extraction into the full manifest.
N=size(data.recordFiles,1); levels={'F0','F1','F2','F3','R','B0'};
primary=find(strcmp(levels,cfg.clinical.featureLevel),1);
for s=1:numel(subjectIndices)
    target=subjectIndices(s);
    if ~strcmp(data.sourceHashes{target},addition.sourceHashes{s})
        error('tbench:InputChanged','Input content/MVC changed during the study. Restart with a fixed input snapshot.');
    end
end
for j=1:numel(candidateIndices)
    k=candidateIndices(j);
    for s=1:numel(subjectIndices), data.recordFiles{subjectIndices(s),k}=addition.recordFiles{s,j}; end
    for l=1:numel(levels)
        incoming=addition.levelViews{l,j};
        if k>size(data.levelViews,2)||isempty(data.levelViews{l,k})
            target=struct('X',nan(N,size(incoming.X,2)),'available',false(N,1), ...
                'summary',repmat(incoming.summary(1),1,N),'featureNames',{incoming.featureNames});
        else, target=data.levelViews{l,k};
        end
        target.X(subjectIndices,:)=incoming.X; target.available(subjectIndices)=incoming.available;
        target.summary(subjectIndices)=incoming.summary; data.levelViews{l,k}=target;
    end
    data.views{k}=data.levelViews{primary,k};
    if isfield(addition,'motionViews')
        incoming=addition.motionViews{j}; old=struct([]);
        if isfield(data,'motionViews')&&k<=numel(data.motionViews)&&~isempty(data.motionViews{k}), old=data.motionViews{k}.rows; end
        % Replace the supplied subjects even when all of their motion is absent.
        if ~isempty(old), old=old(~ismember({old.subject},dataSubjectIds(addition,incoming,subjectIndices))); end
        if isempty(old), old=incoming.rows; elseif ~isempty(incoming.rows), old=[old(:);incoming.rows(:)]; end
        data.motionViews{k}=struct('rows',old);
    end
end
data.cacheHits=data.cacheHits+addition.cacheHits; data.cacheMisses=data.cacheMisses+addition.cacheMisses;
end
function ids=dataSubjectIds(addition,incoming,indices)
ids={};
if ~isempty(incoming.rows), ids={incoming.rows.subject}; end
if ~isempty(addition.levelViews)&&~isempty(addition.levelViews{1,1})
    S=addition.levelViews{1,1}.summary;
    if isfield(S,'subject'), ids=[ids {S.subject}]; end
end
% indices are not subject names; never use them as a guessed identity.
if isempty(ids)&&~isempty(indices)
    error('tbench:MergeMotionIDs','Partial motion merge needs explicit subject IDs.');
end
end
