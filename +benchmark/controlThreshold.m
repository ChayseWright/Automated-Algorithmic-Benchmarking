function [threshold,audit] = controlThreshold(scores,valid,groups,alpha)
%CONTROLTHRESHOLD Strict score > threshold; all tuning controls remain counted.
control=strcmp(groups,'control'); use=control(:)&valid(:)&isfinite(scores(:));
v=sort(scores(use),'descend'); n=numel(v); allowed=floor(alpha*n+1e-12);
threshold=Inf;
if n>0
    if allowed>=n, threshold=-Inf; else, threshold=v(allowed+1); end
end
audit=struct('controlsPlanned',sum(control),'controlsUsed',n,'allowedErrors',allowed, ...
    'empiricalFPR',NaN,'nominalFPR',alpha,'populationGuarantee',false);
if n>0, audit.empiricalFPR=sum(v>threshold)/n; end
end
