function [records,events] = collectRecords(data,candidates,subjects,indices,candidateIndices)
%COLLECTRECORDS Load only selected subjects/candidates; cache identities are retagged.
if isscalar(candidateIndices), candidateIndices=repmat(candidateIndices,size(indices)); end
records=cell(numel(indices),1); events=cell(numel(indices),1);
for j=1:numel(indices)
    s=indices(j); k=candidateIndices(j); file=load(data.recordFiles{s,k},'payload'); p=file.payload;
    records{j}=tag(p.records,candidates(k),subjects(s)); events{j}=tag(p.events,candidates(k),subjects(s));
end
records=join(records); events=join(events);
end
function rows=tag(rows,p,s)
for k=1:numel(rows)
    rows(k).candidateId=p.candidateId; rows(k).candidateHash=p.hash;
    rows(k).pipelineId=p.id; rows(k).preprocessId=p.preprocessId;
    rows(k).group=s.group; rows(k).phenotype=s.phenotype;
    rows(k).severityScore=s.severity_score; rows(k).severityKind=s.severity_kind;
end
end
function result=join(C)
C=C(~cellfun('isempty',C)); result=struct([]);
if ~isempty(C), C=cellfun(@(x)x(:),C,'UniformOutput',false); result=vertcat(C{:}); end
end
