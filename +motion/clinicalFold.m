function [P,saved] = clinicalFold(view,subjects,train,test,cfg,outer)
%CLINICALFOLD Independent motion-only comparator with training-control threshold.
P=struct([]); saved=struct(); if ~cfg.motion.enabled, return; end
S=subjects(train); K=min([cfg.tuning.innerFolds sum(strcmp({S.group},'ET')) sum(strcmp({S.group},'control'))]);
plan=tbench.study.stratifiedFolds(S,K,cfg.tuning.seed+outer); indices=find(train); scores=struct([]);
local=cfg; local.clinical.featureLevel='F1'; local.motion.enabled=false;
for lambda=cfg.tuning.ridgeValues(:)'
    cells=cell(1,K);
    for h=1:K
        fit=indices([plan.fold]~=h); hold=indices([plan.fold]==h);
        model=tbench.benchmark.fitClinical(view.X(fit,:),subjects(fit),view.available(fit),local,lambda);
        cells{h}=tbench.benchmark.predictClinical(model,view.X(hold,:),subjects(hold),view.available(hold),local,.5);
    end
    q=[cells{:}]; [threshold,audit]=tbench.benchmark.controlThreshold([q.score],strcmp({q.status},'ok'),{q.group},cfg.clinical.controlFPR);
    q=tbench.benchmark.applyThreshold(q,threshold,local); metrics=tbench.benchmark.metrics(q);
    violation=max(0,cfg.clinical.minClinicalCoverage-metrics.clinicalCoverage);
    violation=violation+max(0,cfg.clinical.minClinicalCoverage-metrics.controlCoverage);
    violation=violation+max(0,cfg.clinical.minClinicalCoverage-metrics.ETCoverage);
    if isfinite(metrics.controlFPR), violation=violation+max(0,metrics.controlFPR-cfg.clinical.controlFPR); else, violation=violation+1; end
    if metrics.severeSubjects>0, violation=violation+max(0,cfg.clinical.minSevereCapture-metrics.severeCapture);
    elseif cfg.clinical.minSevereCapture>0, violation=violation+1; end
    if strcmp(cfg.clinical.priority,'mild')&&metrics.mildSubjects==0, violation=violation+1; end
    r=struct('eligible',violation<=1e-12,'violation',violation,'metrics',metrics,'ridge',lambda,'threshold',threshold,'thresholdAudit',audit);
    if isempty(scores), scores=r; else, scores(end+1)=r; end %#ok<AGROW>
end
[pick,status]=tbench.study.chooseCandidate(scores,local); selected=scores(pick);
model=tbench.benchmark.fitClinical(view.X(train,:),subjects(train),view.available(train),local,selected.ridge);
P=tbench.benchmark.predictClinical(model,view.X(test,:),subjects(test),view.available(test),local,selected.threshold);
for k=1:numel(P)
    P(k).pipelineId='motion_only_reference'; P(k).level='motion_F1'; P(k).outerFold=outer;
    P(k).selectionStatus=status; P(k).modelRidge=selected.ridge;
end
saved=struct('model',model,'threshold',selected.threshold,'thresholdAudit',selected.thresholdAudit, ...
    'ridge',selected.ridge,'selectionStatus',status,'trainingSubjects',{{subjects(train).id}}, ...
    'targetKind',cfg.motion.clinicalTargetKind,'referencePolicyHash',tbench.motion.policyHash(cfg.motion));
end
