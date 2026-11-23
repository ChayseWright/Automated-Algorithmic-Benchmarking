function [winner,scores,provenance,trace,candidates,data,calibration] = optimize( ...
    subjects,training,initialIndices,candidates,data,calibration,cfg,outerFold,codeHash)
%OPTIMIZE Coverage search then training-only coordinate refinement.
% Candidate pools learned in OTHER folds never enter this fold's search scope.
allowed=initialIndices(:)'; scope='outer_training';
if all(training), scope='all_training_deployment'; end
[winner,scores,provenance]=tbench.study.tune(subjects,training,allowed,candidates,data.views,calibration,cfg,outerFold,data.motionViews);
scores=tagScores(scores,0,scope); provenance=tagProvenance(provenance,0,scope);
[initialWinner,initialStatus]=tbench.study.chooseCandidate(scores,cfg);
winner=scores(initialWinner); winner.selectionStatus=initialStatus;
trace=struct([]); termination='refinement_disabled';
rounds=cfg.tuning.refinementRounds; if ~cfg.tuning.enabled, rounds=0; end
for round=1:rounds
    parent=candidates(winner.candidateIndex); tested={candidates(allowed).hash};
    proposals=struct([]); ranked=scores; parentHashes={};
    for beam=1:cfg.tuning.refinementSeeds
        if isempty(ranked), break; end
        [seedIndex,~]=tbench.study.chooseCandidate(ranked,cfg); seed=candidates(ranked(seedIndex).candidateIndex);
        parentHashes{end+1}=seed.hash; %#ok<AGROW>
        local=cfg; local.tuning.refinementCandidates=max(1,ceil(cfg.tuning.refinementCandidates/cfg.tuning.refinementSeeds));
        seen=tested; if ~isempty(proposals), seen=[seen {proposals.hash}]; end
        neighbors=tbench.study.coordinateNeighbors(seed,local,seen,cfg.tuning.seed+outerFold*1009+round*31+beam);
        proposals=append(proposals,neighbors);
        ranked=ranked(~strcmp({ranked.candidateHash},seed.hash));
    end
    if numel(proposals)>cfg.tuning.refinementCandidates, proposals=proposals(1:cfg.tuning.refinementCandidates); end
    if isempty(proposals), termination='no_unseen_coordinate_neighbor'; break; end
    trialIndices=zeros(1,numel(proposals));
    for j=1:numel(proposals)
        k=find(strcmp({candidates.hash},proposals(j).hash),1);
        if isempty(k), candidates(end+1)=proposals(j); k=numel(candidates); end %#ok<AGROW>
        trialIndices(j)=k;
    end
    % Only outer-training subjects are extracted for learned proposals.
    data=tbench.study.ensureCandidateData(data,subjects,find(training),trialIndices,candidates,cfg,codeHash);
    fresh=proposals;
    if ~isempty(calibration), fresh=fresh(~ismember({fresh.hash},{calibration.candidateHash})); end
    added=tbench.study.calibrationCached(fresh,cfg,codeHash,0);
    calibration=append(calibration,added);
    [~,next,more]=tbench.study.tune(subjects,training,trialIndices,candidates,data.views,calibration,cfg,outerFold,data.motionViews);
    next=tagScores(next,round,scope); more=tagProvenance(more,round,scope);
    scores=append(scores,next); provenance=append(provenance,more);
    previous=winner; [index,status]=tbench.study.chooseCandidate(scores,cfg);
    winner=scores(index); winner.selectionStatus=status; allowed=unique([allowed trialIndices],'stable');
    row=struct('pipelineId',parent.id,'outerFold',outerFold,'scope',scope,'round',round, ...
        'parentCandidateHash',strjoin(parentHashes,' '),'winnerCandidateHash',winner.candidateHash, ...
        'evaluatedSignalConfigurations',numel(allowed),'proposedThisRound',numel(trialIndices), ...
        'winnerChanged',~strcmp(previous.configurationHash,winner.configurationHash), ...
        'trainingSubjects',strjoin({subjects(training).id},' '),'termination','budget_limit');
    trace=append(trace,row); termination='budget_limit';
end
if isempty(trace)
    p=candidates(winner.candidateIndex);
    trace=struct('pipelineId',p.id,'outerFold',outerFold,'scope',scope,'round',0, ...
        'parentCandidateHash',p.hash,'winnerCandidateHash',winner.candidateHash, ...
        'evaluatedSignalConfigurations',numel(allowed),'proposedThisRound',0, ...
        'winnerChanged',false,'trainingSubjects',strjoin({subjects(training).id},' '),'termination',termination);
else, trace(end).termination=termination;
end
winner.searchScope=scope;
end
function rows=tagScores(rows,round,scope)
for k=1:numel(rows), rows(k).searchRound=round; rows(k).selectionScope=scope; end
end
function rows=tagProvenance(rows,round,scope)
for k=1:numel(rows), rows(k).searchRound=round; rows(k).selectionScope=scope; end
end
function a=append(a,b)
if isempty(b), return; end
if isempty(a), a=b(:); else, a=[a(:);b(:)]; end
end
