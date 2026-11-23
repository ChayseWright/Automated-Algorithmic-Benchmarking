function result = inspectSelected(subjects,pipelines,cfg,folder,outerFold)
%INSPECTSELECTED Retain selected signals and compute B5 on outer-test people.
result=struct('prediction',struct([]),'network',struct([]),'frequencySupport',struct([]), ...
    'residualAudit',struct([]),'metrology',struct([]),'temporalPrediction',struct([]),'controlInjection',struct([]),'motionCoherence',struct([]));
algorithmCfg=tbench.motion.inferenceConfig(cfg);
for s=1:numel(subjects)
    motionFiles=cell(numel(cfg.tasks),numel(cfg.trials));
    if cfg.motion.enabled, [motionFiles,~]=tbench.motion.runtimeSubject(cfg,subjects(s).id); end
    try, subject=tbench.loadSubject(subjects(s));
    catch err
        if cfg.failFast, rethrow(err); end
        continue;
    end
    metrologyDone=false;
    for j=1:numel(cfg.tasks)
        task=cfg.tasks{j};
        for r=1:numel(cfg.trials)
            trialName=cfg.trials{r};
            if ~isfield(subject.data,task)||~isfield(subject.data.(task),trialName), continue; end
            try
                trial=tbench.validateTrial(subject.data.(task).(trialName),cfg);
                ctx=tbench.preprocess.reference(trial,subject.mvc,cfg);
            catch err
                if cfg.failFast, rethrow(err); end
                continue;
            end
            outputs=cell(numel(pipelines),15);
            for a=1:numel(pipelines)
                for m=1:15
                    c=tbench.study.channel(ctx,m);
                    try
                        pp=pipelines(a).preprocess(c,algorithmCfg,pipelines(a)); out=pipelines(a).analyze(pp,algorithmCfg,pipelines(a));
                    catch err
                        if cfg.failFast, rethrow(err); end
                        out=tbench.emptyOutput(c); out.status='method_error'; out.message=err.message;
                    end
                    outputs{a,m}=out;
                end
            end
            methods=struct('name',{},'run',{},'params',{});
            for a=1:numel(pipelines)
                methods(a)=struct('name',pipelines(a).id,'run',@tbench.methodWelch,'params',struct());
            end
            local=cfg;
            if cfg.motion.enabled&&cfg.motion.coherence.enabled
                motionFile=load(motionFiles{j,r},'payload'); motionReference=motionFile.payload;
                for a=1:numel(pipelines)
                    packet=[];
                    for m=1:15, packet=tbench.motion.packet(packet,outputs{a,m},m,tbench.study.channel(ctx,m)); end
                    for z=1:numel(motionReference.targets)
                        target=motionReference.targets(z);
                        if ~isempty(cfg.motion.primaryTargets)&&~any(strcmp(cfg.motion.primaryTargets,target.name)), continue; end
                        mPolicy=cfg.motion; mPolicy.sourceSynchronized=cfg.synchronized;
                        detail=tbench.motion.coherence(packet,motionReference,target,mPolicy);
                        for q=1:numel(detail)
                            detail(q).pipelineId=pipelines(a).id; detail(q).candidateId=pipelines(a).candidateId;
                            detail(q).candidateHash=pipelines(a).hash;
                        end
                        detail=tag(detail,subjects(s),task,trialName,outerFold);
                        tbench.study.appendCSV(fullfile(folder,'motion_coherence.csv'),detail);
                        row=struct('pipelineId',pipelines(a).id,'candidateId',pipelines(a).candidateId, ...
                            'target',target.name,'status',detail(1).status,'meanCoherentFraction',finiteMean([detail.coherentFraction]), ...
                            'identifiableInputs',sum([detail.linearModelConditionsMet]),'segments',detail(1).segments);
                        result.motionCoherence=append(result.motionCoherence,tag(row,subjects(s),task,trialName,outerFold));
                    end
                end
            end
            local.prediction.guardSeconds=max(cfg.prediction.guardSeconds,ctx.observerGuardSeconds+1/ctx.fs);
            if cfg.bench.temporal
                T=tbench.study.temporalProbe(trial,subject.mvc,pipelines,local,ctx);
                result.temporalPrediction=append(result.temporalPrediction,tag(T,subjects(s),task,trialName,outerFold));
            end
            if cfg.bench.controlInjection&&strcmp(subjects(s).group,'control')&& ...
                    any(strcmp(cfg.injection.tasks,task))&&any(strcmp(cfg.injection.trials,trialName))
                I=tbench.study.controlInjection(trial,subject.mvc,pipelines,local,outputs,ctx);
                result.controlInjection=append(result.controlInjection,tag(I,subjects(s),task,trialName,outerFold));
            end
            if cfg.execution.computeExplanatory
                if cfg.bench.prediction
                    detail=tbench.predictive(ctx,outputs,methods,local);
                    writeDetail(folder,detail,'prediction',subjects(s),task,trialName,outerFold); result.prediction=append(result.prediction, ...
                        tag(tbench.study.explanatorySummary(detail,'prediction'),subjects(s),task,trialName,outerFold));
                end
                if cfg.bench.network
                    detail=tbench.network(ctx,outputs,methods,local); writeDetail(folder,detail,'network',subjects(s),task,trialName,outerFold);
                    result.network=append(result.network,tag(tbench.study.explanatorySummary(detail,'network'),subjects(s),task,trialName,outerFold));
                    detail=tbench.frequencySupport(ctx,outputs,methods,local); writeDetail(folder,detail,'frequencySupport',subjects(s),task,trialName,outerFold);
                    result.frequencySupport=append(result.frequencySupport,tag(tbench.study.explanatorySummary(detail,'frequencySupport'),subjects(s),task,trialName,outerFold));
                end
                detail=tbench.residualAudit(ctx,outputs,methods,local); writeDetail(folder,detail,'residualAudit',subjects(s),task,trialName,outerFold);
                result.residualAudit=append(result.residualAudit,tag(tbench.study.explanatorySummary(detail,'residualAudit'),subjects(s),task,trialName,outerFold));
            end
            if cfg.execution.saveSelectedSignals
                safeId=regexprep(subjects(s).id,'[^A-Za-z0-9_-]','_');
                basename=[task '_' trialName '.mat'];
                referencePath=fullfile(folder,'references',safeId,basename);
                reference=struct('time',ctx.t,'samplingHz',ctx.fs,'valid',ctx.valid, ...
                    'envelope',ctx.envelope,'energy',ctx.energy,'originalTimeOrigin',ctx.timeOrigin, ...
                    'subject',subjects(s).id,'outerFold',outerFold,'task',task,'trial',trialName, ...
                    'reduction',ctx.reductionAudit,'rawReduction',ctx.rawReductionAudit, ...
                    'runId',cfg.runId,'runNumber',cfg.runNumber);
                if cfg.motion.enabled
                    reference.motionReferenceRelativePath=fullfile('..','..','motion','references',safeId,basename);
                end
                tbench.study.atomicSave(referencePath,reference);
                for a=1:numel(pipelines)
                    destination=fullfile(folder,'selected',sprintf('outer_%02d',outerFold),pipelines(a).id, ...
                        pipelines(a).candidateId,safeId);
                    saved=outputs(a,:);
                    for m=1:numel(saved), saved{m}.time=[]; end % Clock is stored exactly once.
                    payload=struct('pipeline',pipelines(a),'methodNames',{{pipelines(a).id}},'outputs',{saved}, ...
                        'referenceRelativePath',fullfile('..','..','..','..','..','references',safeId,basename));
                    tbench.study.atomicSave(fullfile(destination,basename),payload);
                end
            end
            if cfg.bench.metrology&&~metrologyDone
                M=tbench.study.metrology(trial,subject.mvc,pipelines,local);
                result.metrology=append(result.metrology,tag(M,subjects(s),task,trialName,outerFold)); metrologyDone=true;
            end
        end
    end
end
end
function v=finiteMean(x)
x=x(isfinite(x)); v=NaN; if ~isempty(x), v=mean(x); end
end
function writeDetail(folder,detail,name,subject,task,trial,fold)
detail=tag(detail,subject,task,trial,fold);
tbench.study.appendCSV(fullfile(folder,'explanatory',['B5_' name '_detail.csv']),detail);
end
function a=append(a,b)
if isempty(b), return; end
if isempty(a), a=b(:); else, a=[a(:);b(:)]; end
end
function rows=tag(rows,subject,task,trial,fold)
for k=1:numel(rows)
    rows(k).subject=subject.id; rows(k).group=subject.group; rows(k).phenotype=subject.phenotype;
    rows(k).task=task; rows(k).trial=trial; rows(k).outerFold=fold;
end
end
