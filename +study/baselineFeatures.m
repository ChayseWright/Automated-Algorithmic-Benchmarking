function [X,available,summary,names] = baselineFeatures(records,subjects,cfg)
%BASELINEFEATURES B0 tonic effort/MVC comparator, independent of detection.
names={'tonicEnergy','MVCAmplitude'}; X=nan(numel(subjects),2); available=false(numel(subjects),1);
summary=repmat(struct('subject','','group','','level','B0','coverage',0,'validRecords',0, ...
    'plannedRecords',0,'posturalTasks',0,'kineticTasks',0,'dynamicAvailable',false, ...
    'meanPower',NaN,'upperRecordPower',NaN,'detectionYield',NaN,'rhythmYield4to8',NaN, ...
    'phenotype','','status','indeterminate'),1,numel(subjects));
field='mvcTonicEnergy'; if strcmp(cfg.clinical.powerDomain,'absolute'), field='tonicEnergy'; end
planned=numel(cfg.tasks)*numel(cfg.trials)*15;
for s=1:numel(subjects)
    r=summary(s); r.subject=subjects(s).id; r.group=subjects(s).group; r.phenotype=subjects(s).phenotype; r.plannedRecords=planned;
    F=records(strcmp({records.subject},subjects(s).id));
    if isempty(F), summary(s)=r; continue; end
    usable=[F.inputValidFraction]>0&isfinite([F.(field)])&isfinite([F.mvcReference])&[F.mvcReference]>0;
    r.validRecords=sum(usable); r.coverage=sum([F(usable).inputValidFraction])/planned;
    V=nan(numel(cfg.tasks),2);
    for j=1:numel(cfg.tasks)
        J=F(usable&strcmp({F.task},cfg.tasks{j}));
        if ~isempty(J), V(j,:)=[mean([J.(field)]) median([J.mvcReference])]; end
    end
    post=cellfun(@(x)x(1)=='P',cfg.tasks(:)'); good=all(isfinite(V),2);
    r.posturalTasks=sum(good(post)); r.kineticTasks=sum(good(~post));
    if r.posturalTasks>0&&r.kineticTasks>0
        X(s,:)=.5*(mean(V(post&good',:),1)+mean(V(~post&good',:),1));
        available(s)=r.coverage>=cfg.clinical.minSubjectCoverage;
        if available(s), r.status='ok'; end
    end
    summary(s)=r;
end
end
