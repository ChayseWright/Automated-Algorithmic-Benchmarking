function [roc,severity,withinFold] = curves(P)
%CURVES Descriptive outer-score curves; never used for selecting a winner.
roc=struct([]); severity=struct([]); withinFold=struct([]);
if isempty(P), return; end
ids=unique({P.pipelineId}); levels=unique({P.level});
for a=1:numel(ids)
    for l=1:numel(levels)
        Q=P(strcmp({P.pipelineId},ids{a})&strcmp({P.level},levels{l}));
        if isempty(Q), continue; end
        good=strcmp({Q.status},'ok'); thresholds=[Inf sort(unique([Q(good).score]),'descend') -Inf];
        for z=1:numel(thresholds)
            R=Q;
            for k=1:numel(R), if good(k), R(k).decision=double(R(k).score>thresholds(z)); end; end
            m=tbench.benchmark.metrics(R);
            r=struct('pipelineId',ids{a},'level',levels{l},'threshold',thresholds(z), ...
                'controlFPR',m.controlFPR,'mildCapture',m.mildCapture,'severeCapture',m.severeCapture, ...
                'ETCapture',m.ETCapture,'controls',m.controls,'ET',m.ET,'coverage',m.clinicalCoverage);
            roc=append(roc,r);
        end
        scores=sort(unique([Q(strcmp({Q.group},'ET')).trueSeverityScore])); scores=scores(isfinite(scores));
        for z=1:numel(scores)
            low=strcmp({Q.group},'ET')&[Q.trueSeverityScore]<=scores(z);
            high=strcmp({Q.group},'ET')&[Q.trueSeverityScore]>=scores(z);
            r=struct('pipelineId',ids{a},'level',levels{l},'severityCut',scores(z), ...
                'lowerCount',sum(low),'lowerCaptured',sum(low&[Q.decision]==1), ...
                'higherCount',sum(high),'higherCaptured',sum(high&[Q.decision]==1));
            severity=append(severity,r);
        end
        folds=unique([Q.outerFold]);
        for z=1:numel(folds)
            m=tbench.benchmark.metrics(Q([Q.outerFold]==folds(z)));
            r=struct('pipelineId',ids{a},'level',levels{l},'outerFold',folds(z), ...
                'AUCavailable',m.AUCavailable,'mildAUCavailable',m.mildAUCavailable, ...
                'severeAUCavailable',m.severeAUCavailable,'controls',m.controls,'ET',m.ET);
            withinFold=append(withinFold,r);
        end
    end
end
end
function a=append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
