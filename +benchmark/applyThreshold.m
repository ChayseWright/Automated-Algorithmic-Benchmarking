function P = applyThreshold(P,threshold,cfg)
%APPLYTHRESHOLD Explicit control operating point; indeterminate stays distinct.
for k=1:numel(P)
    P(k).threshold=threshold;
    if ~strcmp(P(k).status,'ok')||~isfinite(threshold)
        P(k).decision=NaN; P(k).severityDecision=NaN; P(k).status='indeterminate'; continue;
    end
    P(k).decision=double(P(k).score>threshold); P(k).severityDecision=0;
    if P(k).decision==1
        p=P(k).score; joint=[1-p p*[P(k).ordinalP1 P(k).ordinalP2 P(k).ordinalP3]];
        costs=joint*cfg.clinical.costMatrix; [~,j]=min(costs(2:4)); P(k).severityDecision=j;
    end
end
end
