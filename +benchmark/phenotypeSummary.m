function rows = phenotypeSummary(P)
%PHENOTYPESUMMARY B8 subsets of an inclusively trained clinical procedure.
rows=struct([]); if isempty(P), return; end
ids=unique({P.pipelineId}); levels=unique({P.level});
for a=1:numel(ids)
    for l=1:numel(levels)
        Q=P(strcmp({P.pipelineId},ids{a})&strcmp({P.level},levels{l}));
        if isempty(Q), continue; end
        names=unique({Q(strcmp({Q.group},'ET')).phenotype});
        for j=1:numel(names)
            use=strcmp({Q.group},'control')|strcmp({Q.phenotype},names{j});
            r=tbench.benchmark.metrics(Q(use)); r.pipelineId=ids{a}; r.level=levels{l};
            r.phenotype=names{j}; r.trainingScope='inclusive_training_subset_evaluation';
            rows=append(rows,r);
        end
    end
end
end
function a=append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
