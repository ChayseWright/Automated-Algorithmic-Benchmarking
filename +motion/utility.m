function s = utility(rows,cfg)
%UTILITY Macro subject, task-type and target weighting; no cross-unit sums.
s=struct('gain',NaN,'loss',NaN,'baselineLoss',NaN,'subjects',0,'inputCoverage',0);
if isempty(rows), return; end
ids=unique({rows.subject}); gain=[]; losses=[]; baselines=[]; coverage=[];
for k=1:numel(ids)
    Q=rows(strcmp({rows.subject},ids{k})); targets=unique({Q.target}); targetGain=[]; targetLoss=[]; targetBase=[]; targetCover=[];
    for j=1:numel(targets)
        P=Q(strcmp({Q.target},targets{j})); typeLoss=[]; typeBase=[]; typeCover=[];
        for type='PK'
            tasks=cfg.tasks(cellfun(@(x)x(1)==type,cfg.tasks)); a=[]; b=[]; c=[];
            for z=1:numel(tasks)
                R=P(strcmp({P.task},tasks{z})); if isempty(R), continue; end
                a(end+1)=mean([R.transformedLoss]); b(end+1)=mean([R.baselineTransformedLoss]); c(end+1)=mean([R.inputAvailable]); %#ok<AGROW>
            end
            if ~isempty(a), typeLoss(end+1)=mean(a); typeBase(end+1)=mean(b); typeCover(end+1)=mean(c); end %#ok<AGROW>
        end
        if numel(typeLoss)==2
            l=mean(typeLoss); b=mean(typeBase);
            if b>cfg.motion.minimumBaselineLoss, targetGain(end+1)=(b-l)/b; end %#ok<AGROW>
            targetLoss(end+1)=l; targetBase(end+1)=b; targetCover(end+1)=mean(typeCover); %#ok<AGROW>
        end
    end
    if ~isempty(targetGain), gain(end+1)=mean(targetGain); end %#ok<AGROW>
    if ~isempty(targetLoss), losses(end+1)=mean(targetLoss); baselines(end+1)=mean(targetBase); coverage(end+1)=mean(targetCover); end %#ok<AGROW>
end
s.subjects=numel(losses);
if ~isempty(gain), s.gain=mean(gain); end
if ~isempty(losses), s.loss=mean(losses); s.baselineLoss=mean(baselines); s.inputCoverage=mean(coverage); end
end
