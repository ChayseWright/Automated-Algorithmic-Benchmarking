function rows = commonEligible(P,cfg)
%COMMONELIGIBLE Secondary shared-output cohort, beside all-enrolled results.
rows=struct([]); if isempty(P), return; end
Q=P(strcmp({P.level},cfg.clinical.featureLevel)); ids=unique({Q.pipelineId});
if isempty(ids), return; end
common=unique({Q.subject});
for a=1:numel(ids)
    A=Q(strcmp({Q.pipelineId},ids{a})&strcmp({Q.status},'ok'));
    common=intersect(common,{A.subject});
end
for a=1:numel(ids)
    A=Q(strcmp({Q.pipelineId},ids{a})&ismember({Q.subject},common));
    r=tbench.benchmark.metrics(A); r.pipelineId=ids{a}; r.level=cfg.clinical.featureLevel;
    r.commonSubjects=numel(common); r.enrolledSubjects=numel(unique({Q.subject}));
    if isempty(rows), rows=r; else, rows(end+1)=r; end %#ok<AGROW>
end
end
