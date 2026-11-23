function [spaces,stability] = parameterTables(study)
%PARAMETERTABLES Declared search and outer-fold parameter stability, separately.
spaces=struct([]); stability=struct([]); cfg=study.config;
for a=1:numel(cfg.pipelines)
    p=cfg.pipelines(a); [pp,ap,defaults,description]=tbench.study.parameterSpace(p,cfg);
    keys=fieldnames(defaults);
    for j=1:numel(keys), if ~isfield(p.analysisParams,keys{j}), p.analysisParams.(keys{j})=defaults.(keys{j}); end; end
    p=tbench.study.finalizeCandidate(p,cfg);
    for module={'preprocessParams','analysisParams','clinical'}
        name=module{1};
        if strcmp(name,'preprocessParams'), space=pp;
        elseif strcmp(name,'analysisParams'), space=ap;
        else, space=struct('ridge',cfg.tuning.ridgeValues);
        end
        keys=fieldnames(space);
        for j=1:numel(keys)
            key=keys{j}; value=[];
            if strcmp(name,'clinical'), value=cfg.clinical.ridge;
            elseif isfield(p.(name),key), value=p.(name).(key);
            end
            row=struct('pipelineId',p.id,'analysisId',p.analysisId,'inputKind',p.inputKind, ...
                'module',name,'parameter',key,'defaultJSON',jsonencode(value), ...
                'choicesJSON',jsonencode(space.(key)),'searchStrategy',description.strategy, ...
                'optimality',description.optimality,'transformPolicy',description.transformPolicy);
            spaces=append(spaces,row);
        end
    end
    selections=study.selection(strcmp({study.selection.pipelineId},p.id));
    if isempty(selections), continue; end
    chosen=study.candidateCatalog([selections.candidateIndex]);
    for module={'preprocessParams','analysisParams','clinical'}
        name=module{1};
        if strcmp(name,'clinical'), keys={'ridge'}; else, keys=fieldnames(chosen(1).(name)); end
        for j=1:numel(keys)
            key=keys{j}; values=cell(1,numel(chosen));
            for k=1:numel(chosen)
                if strcmp(name,'clinical'), value=selections(k).ridge; else, value=chosen(k).(name).(key); end
                values{k}=jsonencode(value);
            end
            uniqueValues=unique(values,'stable');
            for k=1:numel(uniqueValues)
                count=sum(strcmp(values,uniqueValues{k}));
                row=struct('pipelineId',p.id,'module',name,'parameter',key, ...
                    'valueJSON',uniqueValues{k},'selectedOuterFolds',count, ...
                    'outerFolds',numel(chosen),'selectionFraction',count/numel(chosen), ...
                    'scope','descriptive_overlapping_training_sets_not_independent_replicates');
                stability=append(stability,row);
            end
        end
    end
end
end
function a=append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
