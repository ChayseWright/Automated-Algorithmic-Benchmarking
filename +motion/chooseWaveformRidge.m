function policy = chooseWaveformRidge(data,index,subjects,training,cfg)
%CHOOSEWAVEFORMRIDGE Training-person repetition validation; no outer-test targets.
policy=struct('ridge',cfg.motion.readoutRidges(1),'gain',NaN,'subjects',0, ...
    'scope','outer_training_people_person_calibrated_repetition_validation', ...
    'trainingSubjects',strjoin({subjects(training).id},' '));
scores=nan(sum(training),numel(cfg.motion.readoutRidges)); n=0;
for s=find(training(:))'
    n=n+1; loaded=load(data.recordFiles{s,index},'payload');
    if ~isfield(loaded.payload,'motion'), continue; end
    rows=loaded.payload.motion.waveformPrediction;
    for k=1:numel(cfg.motion.readoutRidges)
        if isempty(rows), continue; end
        Q=rows([rows.ridge]==cfg.motion.readoutRidges(k)&isfinite([rows.gainOverNull]));
        values=[];
        for type='PK'
            tasks=cfg.tasks(cellfun(@(x)x(1)==type,cfg.tasks)); v=[];
            for j=1:numel(tasks)
                R=Q(strcmp({Q.task},tasks{j})); if ~isempty(R), v(end+1)=mean([R.gainOverNull]); end %#ok<AGROW>
            end
            if ~isempty(v), values(end+1)=mean(v); end %#ok<AGROW>
        end
        if numel(values)==2, scores(n,k)=mean(values); end
    end
end
means=nan(1,size(scores,2)); counts=zeros(size(means));
for k=1:numel(means), q=scores(:,k); q=q(isfinite(q)); counts(k)=numel(q); if ~isempty(q), means(k)=mean(q); end; end
means(~isfinite(means))=-Inf; [best,pick]=max(means);
if isfinite(best), policy.ridge=cfg.motion.readoutRidges(pick); policy.gain=best; policy.subjects=counts(pick); end
end
