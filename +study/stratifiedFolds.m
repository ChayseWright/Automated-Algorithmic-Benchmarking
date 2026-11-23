function plan = stratifiedFolds(subjects,K,seed)
%STRATIFIEDFOLDS Group by participant; balance controls and ordered ET severity.
state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
rng(seed,'twister'); N=numel(subjects); assignment=zeros(1,N);
for group={'control','ET'}
    ix=find(strcmp({subjects.group},group{1}));
    if isempty(ix), continue; end
    if strcmp(group{1},'ET')
        severity=[subjects(ix).severity_score]; jitter=rand(size(severity))*1e-10;
        severity(~isfinite(severity))=Inf; [~,q]=sort(severity+jitter); ix=ix(q);
        % Consecutive severity blocks distribute across all folds.
        for start=1:K:numel(ix)
            block=ix(start:min(numel(ix),start+K-1)); order=randperm(K);
            assignment(block)=order(1:numel(block));
        end
    else
        ix=ix(randperm(numel(ix))); assignment(ix)=mod(0:numel(ix)-1,K)+1;
    end
end
plan=repmat(struct('subject','','fold',0,'group',''),1,N);
for k=1:N, plan(k)=struct('subject',subjects(k).id,'fold',assignment(k),'group',subjects(k).group); end
end
