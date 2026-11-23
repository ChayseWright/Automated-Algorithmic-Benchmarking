function prediction = predictClinical(model,X,subjects,available,cfg,threshold)
%PREDICTCLINICAL Case flag and severity decision use explicit separate policies.
for k=1:numel(model.featureTransform.columns)
    j=model.featureTransform.columns(k); X(:,j)=asinh(X(:,j)/model.featureTransform.scales(k));
end
p=tbench.benchmark.predictLogistic(model.diagnosis,X);
fallback=~available(:)|~all(isfinite(X),2); p(fallback)=model.trainingPrior;
n=size(X,1); conditional=ones(n,3)/3;
if ~isempty(model.cumulative{1})&&~isempty(model.cumulative{2})
    c1=tbench.benchmark.predictLogistic(model.cumulative{1},X);
    c2=tbench.benchmark.predictLogistic(model.cumulative{2},X);
    c2=max(c1,c2); conditional=[c1 c2-c1 1-c2];
end
joint=[1-p bsxfun(@times,p,conditional)]; numeric=nan(n,1);
if ~isempty(model.numericSeverity)
    s=model.numericSeverity; Z=X;
    for k=1:size(Z,2), Z(~isfinite(Z(:,k)),k)=model.diagnosis.fill(k); end
    B=[ones(n,1) bsxfun(@rdivide,bsxfun(@minus,Z,s.mu),s.sd)];
    numeric=max(s.minimum,min(s.maximum,B*s.beta));
end
prediction=repmat(struct('subject','','group','','score',NaN,'decision',NaN, ...
    'severityPrediction',NaN,'severityDecision',NaN,'trueSeverityGroup',0, ...
    'trueSeverityScore',NaN,'severityKind','','severityDefinition','', ...
    'ordinalP1',NaN,'ordinalP2',NaN,'ordinalP3',NaN,'status','indeterminate', ...
    'threshold',threshold,'trainingPrior',model.trainingPrior,'phenotype',''),1,n);
% Applying training cuts to test labels is allowed for evaluation, not fitting.
local=cfg; local.clinical.severityCuts=model.severityCuts;
if isempty(model.severityCuts)&&~strcmp(cfg.clinical.severityMode,'provided')
    trueGroup=zeros(1,n);
else, trueGroup=tbench.benchmark.severityGroups(subjects,local,false(1,n));
end
cost=joint*cfg.clinical.costMatrix;
for k=1:n
    r=prediction(k); r.subject=subjects(k).id; r.group=subjects(k).group; r.score=p(k);
    r.trueSeverityScore=subjects(k).severity_score; r.trueSeverityGroup=trueGroup(k);
    r.severityKind=subjects(k).severity_kind; r.severityDefinition=model.severityDescription;
    r.severityPrediction=numeric(k); r.ordinalP1=conditional(k,1); r.ordinalP2=conditional(k,2); r.ordinalP3=conditional(k,3);
    r.phenotype=subjects(k).phenotype;
    if ~fallback(k)&&isfinite(threshold)&&strcmp(model.status,'ok')
        r.status='ok'; r.decision=double(p(k)>threshold); r.severityDecision=0;
        if r.decision==1
            [~,j]=min(cost(k,2:4)); r.severityDecision=j;
        end
    end
    prediction(k)=r;
end
end
