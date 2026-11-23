function [P,summary,model,thresholdAudit] = evaluateLevel(records,subjects,train,test,cfg,level,lambda,outerFold)
%EVALUATELEVEL Training-calibrated compression/readout without test-label tuning.
if nargin<8, outerFold=0; end
local=cfg; local.clinical.featureLevel=level;
if isstruct(records)&&isscalar(records)&&isfield(records,'X')
    X=records.X; available=records.available; summary=records.summary;
else
    [X,available,summary]=tbench.study.subjectFeatures(records,subjects,local,level);
end
ts=subjects(train); nc=sum(strcmp({ts.group},'control')); ne=sum(strcmp({ts.group},'ET'));
K=min([cfg.tuning.innerFolds nc ne]); cells=cell(1,K);
plan=tbench.study.stratifiedFolds(ts,K,cfg.tuning.seed+outerFold);
indices=find(train);
for h=1:K
    fit=indices([plan.fold]~=h); hold=indices([plan.fold]==h);
    inner=tbench.benchmark.fitClinical(X(fit,:),subjects(fit),available(fit),local,lambda);
    cells{h}=tbench.benchmark.predictClinical(inner,X(hold,:),subjects(hold),available(hold),local,.5);
end
Q=[cells{:}];
[threshold,thresholdAudit]=tbench.benchmark.controlThreshold([Q.score],strcmp({Q.status},'ok'),{Q.group},cfg.clinical.controlFPR);
model=tbench.benchmark.fitClinical(X(train,:),subjects(train),available(train),local,lambda);
P=tbench.benchmark.predictClinical(model,X(test,:),subjects(test),available(test),local,threshold);
end
