function result = evaluateFold(view,subjects,train,test,cfg,outer,pipeline,winner)
%EVALUATEFOLD Population motion prediction/retention on untouched subjects.
result=struct('predictions',struct([]),'policies',struct([]),'models',struct([]));
if ~cfg.motion.enabled||isempty(view.rows), return; end
S=subjects(train); K=min([cfg.tuning.innerFolds sum(strcmp({S.group},'ET')) sum(strcmp({S.group},'control'))]);
plan=tbench.study.stratifiedFolds(S,K,cfg.tuning.seed+outer);
levels={'F0','F1','R'};
if pipeline.expectedDynamic, levels={'F0','F1','F2','F3','R'}; end
[grades,~,description]=tbench.benchmark.severityGroups(subjects,cfg,train);
richRidge=winner.motion.ridge;
for j=1:numel(levels)
    level=levels{j}; policy=tbench.motion.trainingScore(view,subjects,train,plan,cfg,level);
    if strcmp(level,'R'), richRidge=policy.ridge; end
    result.policies=append(result.policies,policyRow(policy,level,pipeline,outer,subjects,train));
    [rows,models]=tbench.motion.lossRows(view,{subjects(train).id},{subjects(test).id},cfg,level,policy.ridge);
    result.models=append(result.models,models);
    rows=tag(rows,subjects,grades,description,pipeline,outer); result.predictions=append(result.predictions,rows);
end
for j=1:numel(levels)
    if strcmp(levels{j},'R'), continue; end
    level=['fixed_' levels{j}];
    [rows,models]=tbench.motion.lossRows(view,{subjects(train).id},{subjects(test).id},cfg,level,richRidge);
    result.models=append(result.models,models);
    rows=tag(rows,subjects,grades,description,pipeline,outer); result.predictions=append(result.predictions,rows);
end
end
function p=policyRow(s,level,pipeline,outer,subjects,train)
p=s; p.level=level; p.pipelineId=pipeline.id; p.candidateId=pipeline.candidateId; p.candidateHash=pipeline.hash;
p.outerFold=outer; p.trainingSubjects=strjoin({subjects(train).id},' ');
end
function rows=tag(rows,subjects,grades,description,p,outer)
for k=1:numel(rows)
    s=find(strcmp({subjects.id},rows(k).subject),1);
    rows(k).group=subjects(s).group; rows(k).phenotype=subjects(s).phenotype;
    rows(k).severityScore=subjects(s).severity_score; rows(k).severityGroup=grades(s); rows(k).severityDefinition=description;
    rows(k).pipelineId=p.id; rows(k).candidateId=p.candidateId; rows(k).candidateHash=p.hash; rows(k).outerFold=outer;
end
end
function a=append(a,b)
if isempty(b), return; end
if isempty(a), a=b(:); else, a=[a(:);b(:)]; end
end
