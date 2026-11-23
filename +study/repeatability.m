function result = repeatability(features,cfg)
%REPEATABILITY Do not let case/control heterogeneity inflate a pooled ICC.
result=struct('repeatability',struct([]),'ICC',struct([]),'topography',struct([]));
if isempty(features), return; end
pipelines=unique({features.pipelineId}); groups=unique({features.group});
for a=1:numel(pipelines)
    for g=1:numel(groups)
        F=features(strcmp({features.pipelineId},pipelines{a})&strcmp({features.group},groups{g}));
        if isempty(F), continue; end
        [R,I,T]=tbench.repeatability(F,cfg);
        result.repeatability=append(result.repeatability,tag(R,pipelines{a},groups{g},'absolute'));
        result.ICC=append(result.ICC,tag(I,pipelines{a},groups{g},'absolute'));
        result.topography=append(result.topography,tag(T,pipelines{a},groups{g},'mvc_composition'));
        normalized=F;
        for k=1:numel(normalized)
            normalized(k).powerMean=normalized(k).mvcPowerMean;
            if ~isfinite(normalized(k).powerMean)&&strcmp(normalized(k).status,'ok'), normalized(k).status='mvc_unavailable'; end
        end
        [R,I,~]=tbench.repeatability(normalized,cfg);
        result.repeatability=append(result.repeatability,tag(R,pipelines{a},groups{g},'mvc'));
        result.ICC=append(result.ICC,tag(I,pipelines{a},groups{g},'mvc'));
    end
end
end
function rows=tag(rows,pipeline,group,domain)
for k=1:numel(rows), rows(k).pipelineId=pipeline; rows(k).group=group; rows(k).powerDomain=domain; end
end
function a=append(a,b)
if isempty(b), return; end
if isempty(a), a=b(:); else, a=[a(:);b(:)]; end
end
