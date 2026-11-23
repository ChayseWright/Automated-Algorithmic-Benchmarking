function candidates = candidates(pipelines,cfg)
%CANDIDATES Exhaust small spaces; otherwise cover parameters, then fill space.
state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
candidates=struct([]);
for k=1:numel(pipelines)
    base=pipelines(k); base.implementationHash=tbench.study.moduleHash(base);
    [pp,ap,defaults,description]=tbench.study.parameterSpace(base,cfg);
    names=fieldnames(defaults);
    for j=1:numel(names), if ~isfield(base.analysisParams,names{j}), base.analysisParams.(names{j})=defaults.(names{j}); end; end
    base.searchDescription=description;
    base=tbench.study.finalizeCandidate(base,cfg); family=base;
    if ~cfg.tuning.enabled, candidates=join(candidates,family); continue; end
    [modules,fields,values]=flatten(pp,ap); total=prod(cellfun(@numel,values));
    budget=cfg.tuning.maxCandidates;
    rng(cfg.tuning.seed+sum(double(base.analysisId).*(1:numel(base.analysisId))),'twister');
    if strcmp(cfg.tuning.searchMode,'grid')||total<=budget
        if total>cfg.tuning.maxGridCandidates
            error('tbench:GridSize','%s requires %.0f configurations; raise maxGridCandidates explicitly or use coverage search.',base.analysisId,total);
        end
        for index=0:total-1
            p=base; number=index;
            for j=1:numel(fields)
                q=mod(number,numel(values{j}))+1; number=floor(number/numel(values{j}));
                p.(modules{j}).(fields{j})=values{j}{q};
            end
            p=tbench.study.finalizeCandidate(p,cfg); family=uniqueAppend(family,p);
        end
    else
        % Main-effect anchors guarantee each declared choice is represented.
        for j=1:numel(fields)
            for q=1:numel(values{j})
                p=base; p.(modules{j}).(fields{j})=values{j}{q};
                p=tbench.study.finalizeCandidate(p,cfg); family=uniqueAppend(family,p);
            end
        end
        if numel(family)>budget
            if cfg.tuning.requireParameterCoverage
                error('tbench:SearchBudget','%s needs at least %d candidates to cover all declared parameter values; requested %d.',base.analysisId,numel(family),budget);
            end
            family=family(1:budget);
        end
        attempts=0;
        while numel(family)<budget&&attempts<100*budget
            attempts=attempts+1; p=base;
            for j=1:numel(fields), p.(modules{j}).(fields{j})=values{j}{randi(numel(values{j}))}; end
            p=tbench.study.finalizeCandidate(p,cfg); family=uniqueAppend(family,p);
        end
    end
    candidates=join(candidates,family);
end
end
function [modules,fields,values]=flatten(pp,ap)
modules={}; fields={}; values={};
for module={'preprocessParams','analysisParams'}
    space=pp; if strcmp(module{1},'analysisParams'), space=ap; end
    names=fieldnames(space);
    for j=1:numel(names)
        v=space.(names{j}); if ~iscell(v), v=num2cell(v); end
        v=v(:)'; if isempty(v), error('tbench:SearchSpace','Empty parameter options.'); end
        modules{end+1}=module{1}; fields{end+1}=names{j}; values{end+1}=v; %#ok<AGROW>
    end
end
end
function a=uniqueAppend(a,b)
if isempty(a)||~any(strcmp({a.hash},b.hash)), a=join(a,b); end
end
function a=join(a,b)
if isempty(a), a=b; else, a=[a(:)' b(:)']; end
end
