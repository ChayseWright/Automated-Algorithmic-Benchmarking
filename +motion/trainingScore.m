function score = trainingScore(view,subjects,training,plan,cfg,level)
%TRAININGSCORE Choose mechanical readout ridge on inner held-out SUBJECTS.
if nargin<6, level='F1'; end
score=struct('available',false,'gain',NaN,'ridge',cfg.motion.readoutRidges(1), ...
    'loss',NaN,'baselineLoss',NaN,'subjects',0,'inputCoverage',0,'scope','inner_subject_validation');
if ~cfg.motion.enabled||isempty(view)||isempty(view.rows), return; end
S=subjects(training); candidates=cfg.motion.readoutRidges; summaries=cell(size(candidates));
for k=1:numel(candidates)
    rows=struct([]);
    for h=1:max([plan.fold])
        fit={S([plan.fold]~=h).id}; hold={S([plan.fold]==h).id};
        q=tbench.motion.lossRows(view,fit,hold,cfg,level,candidates(k));
        if isempty(rows), rows=q; elseif ~isempty(q), rows=[rows(:);q(:)]; end
    end
    summaries{k}=tbench.motion.utility(rows,cfg);
end
gain=cellfun(@(x)x.gain,summaries); gain(~isfinite(gain))=-Inf; [best,index]=max(gain);
if ~isfinite(best), return; end
s=summaries{index}; names=fieldnames(s);
for k=1:numel(names), score.(names{k})=s.(names{k}); end
score.available=true; score.ridge=candidates(index);
end
