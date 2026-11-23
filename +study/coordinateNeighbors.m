function neighbors = coordinateNeighbors(parent,cfg,alreadyTested,seed)
%COORDINATENEIGHBORS One coordinate change at the current training-only winner.
state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
rng(seed,'twister'); [pp,ap]=tbench.study.parameterSpace(parent,cfg);
neighbors=struct([]); modules={}; fields={}; values={};
for module={'preprocessParams','analysisParams'}
    space=pp; if strcmp(module{1},'analysisParams'), space=ap; end
    names=fieldnames(space);
    for j=1:numel(names)
        v=space.(names{j}); if ~iscell(v), v=num2cell(v); end
        modules{end+1}=module{1}; fields{end+1}=names{j}; values{end+1}=v(:)'; %#ok<AGROW>
    end
end
% Randomized dimension order is deterministic and does not inspect test labels.
order=randperm(numel(fields)); queues=cell(size(order)); longest=0;
for j=1:numel(order)
    d=order(j); options=values{d}; options=options(randperm(numel(options))); queue=struct([]);
    for v=1:numel(options)
        p=parent; p.(modules{d}).(fields{d})=options{v}; p=tbench.study.finalizeCandidate(p,cfg);
        if any(strcmp(alreadyTested,p.hash)), continue; end
        if ~isempty(queue)&&any(strcmp({queue.hash},p.hash)), continue; end
        if isempty(queue), queue=p; else, queue(end+1)=p; end %#ok<AGROW>
    end
    queues{j}=queue; longest=max(longest,numel(queue));
end
% Round robin across dimensions prevents one long coordinate dominating a budget.
for depth=1:longest
    for j=1:numel(queues)
        Q=queues{j}; if depth>numel(Q), continue; end
        p=Q(depth);
        if ~isempty(neighbors)&&any(strcmp({neighbors.hash},p.hash)), continue; end
        if isempty(neighbors), neighbors=p; else, neighbors(end+1)=p; end %#ok<AGROW>
        if numel(neighbors)>=cfg.tuning.refinementCandidates, return; end
    end
end
end
