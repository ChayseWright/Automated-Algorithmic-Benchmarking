function data = extract(subjects,candidates,cfg,codeHash)
%EXTRACT One raw subject load; shared preprocessing; validated feature caches.
N=numel(subjects); K=numel(candidates); perSubject=numel(cfg.tasks)*numel(cfg.trials)*15;
algorithmCfg=tbench.motion.inferenceConfig(cfg);
motionEnabled=isfield(cfg,'motion')&&cfg.motion.enabled;
batch=tbench.option(cfg.execution,'candidateBatchSize',32);
if K>batch
    % Bound retained per-subject feature structs. This trades extra file reads
    % for bounded RAM; every original sample still enters each selected recipe.
    for first=1:batch:K
        ix=first:min(K,first+batch-1);
        if cfg.verbose, fprintf('Candidate batch %d--%d / %d\n',ix(1),ix(end),K); end
        part=tbench.study.extract(subjects,candidates(ix),cfg,codeHash);
        if first==1, data=part;
        else, data=tbench.study.mergeCandidateData(data,part,1:N,ix,cfg);
        end
    end
    return;
end
inventories=cell(N,1);
data=struct('recordFiles',{cell(N,K)},'views',{cell(1,K)},'levelViews',{cell(6,K)},'inventory',struct([]), ...
    'cacheHits',0,'cacheMisses',0,'sourceHashes',{{}},'reduction',struct([]), ...
    'motionViews',{repmat({struct('rows',struct([]))},1,K)});
sourceHashes=cell(N,1); processing=tbench.study.processingConfig(cfg);
ppKeys=cell(1,K); detectorKeys=cell(1,K);
for k=1:K
    ppKeys{k}=tbench.study.preparationKey(candidates(k)); p=candidates(k).preprocessParams;
    detectorKeys{k}=sprintf('%.17g_%.17g_%.17g',p.rawHighpassHz,p.envelopePassHz,p.envelopeStopHz);
end
[uniqueDetector,~,detectorIndex]=unique(detectorKeys);
[uniquePP,~,ppIndex]=unique(ppKeys); %#ok<ASGLU>
levels={'F0','F1','F2','F3','R','B0'};
for k=1:K
    for l=1:numel(levels)
        data.levelViews{l,k}=struct('X',[],'available',false(N,1),'summary',struct([]),'featureNames',{{}});
    end
end
cacheRoot=fullfile(cfg.cacheDir,'features'); if ~exist(cacheRoot,'dir'), mkdir(cacheRoot); end
for s=1:N
    if cfg.verbose, fprintf('Extract %s (%d/%d), %d candidates\n',subjects(s).id,s,N,K); end
    sourceError='';
    motionFiles=cell(numel(cfg.tasks),numel(cfg.trials)); motionKey='motion_disabled';
    if motionEnabled, [motionFiles,motionKey]=tbench.motion.runtimeSubject(cfg,subjects(s).id); end
    try, sourceHashes{s}=tbench.study.sourceHash(subjects(s));
    catch err
        if cfg.failFast, rethrow(err); end
        sourceError=err.message; sourceHashes{s}=['unavailable_' tbench.study.hash(struct('id',subjects(s).id,'error',err.message))];
    end
    missing=true(1,K); keys=cell(1,K); paths=cell(1,K); redLoaded=false;
    for k=1:K
        ap=candidates(k).analysisParams; if isfield(ap,'ridge'), ap=rmfield(ap,'ridge'); end
        signalSignature=struct('source',sourceHashes{s},'subject',subjects(s).id,'pipelineId',candidates(k).id,'code',codeHash, ...
            'config',processing,'implementation',candidates(k).implementationHash, ...
            'preprocess',candidates(k).preprocessId,'preprocessOwner',candidates(k).preprocessAnalysisId, ...
            'pp',candidates(k).preprocessParams, ...
            'analysis',candidates(k).analysisId,'ap',ap);
        keys{k}=tbench.study.hash(signalSignature); paths{k}=fullfile(cacheRoot,[keys{k} '.mat']);
        if cfg.execution.cache&&cfg.execution.resumeCache&&exist(paths{k},'file')&&isempty(sourceError)
            try
                file=load(paths{k},'payload'); p=file.payload;
                motionOK=~motionEnabled||(isfield(p,'motionKey')&&strcmp(p.motionKey,motionKey)&&isfield(p,'motion'));
                if strcmp(p.key,keys{k})&&strcmp(p.sourceHash,sourceHashes{s})&&numel(p.records)==perSubject&&motionOK
                    data.recordFiles{s,k}=paths{k}; inventories{s}=p.inventory;
                    data=addViews(data,p.records,subjects(s),cfg,levels,s,k,N);
                    if motionEnabled, data.motionViews{k}.rows=append(data.motionViews{k}.rows,p.motion.featureRows); end
                    if ~redLoaded&&isfield(p,'reduction')
                        data.reduction=append(data.reduction,p.reduction); redLoaded=true;
                    end
                    missing(k)=false; data.cacheHits=data.cacheHits+1;
                end
            catch
                % An incomplete/corrupt cache is never accepted as a completed result.
            end
        end
    end
    if ~any(missing), continue; end
    data.cacheMisses=data.cacheMisses+sum(missing);
    try
        if ~isempty(sourceError), error('tbench:Source','%s',sourceError); end
        subject=tbench.loadSubject(subjects(s)); sourceError='';
    catch err
        if cfg.failFast, rethrow(err); end
        subject=struct('data',struct(),'mvc',nan(1,15)); sourceError=err.message;
    end
    rows=cell(1,K); events=cell(1,K); inv=struct([]); reduction=struct([]); index=0;
    motionResults=repmat({struct('agreement',struct([]),'waveformPrediction',struct([]),'featureRows',struct([]))},1,K);
    for j=1:numel(cfg.tasks)
        task=cfg.tasks{j};
        taskPackets=cell(K,numel(cfg.trials)); taskShared=cell(1,numel(cfg.trials));
        for r=1:numel(cfg.trials)
            trialName=cfg.trials{r}; validTrial=true; message=''; rawFs=NaN; duration=NaN; ctx=[];
            try
                if ~isempty(sourceError), error('tbench:Source','%s',sourceError); end
                if ~isfield(subject.data,task)||~isfield(subject.data.(task),trialName)
                    error('tbench:Missing','Missing %s.%s.',task,trialName);
                end
                trial=tbench.validateTrial(subject.data.(task).(trialName),cfg); rawFs=trial.fs; duration=trial.duration;
                ctx=tbench.preprocess.reference(trial,subject.mvc,cfg);
            catch err
                if cfg.failFast, rethrow(err); end
                validTrial=false; message=err.message;
            end
            status='ok'; if ~validTrial, status='invalid'; end
            invRow=struct('subject',subjects(s).id,'task',task,'trial',trialName,'status',status, ...
                'message',message,'rawSamplingHz',rawFs,'durationSeconds',duration);
            inv=append(inv,invRow);
            if validTrial
                red=ctx.reductionAudit; red.subject=subjects(s).id; red.task=task; red.trial=trialName;
                red.carrierOutputFs=ctx.reducedRawFs; reduction=append(reduction,red);
            end
            if motionEnabled
                file=load(motionFiles{j,r},'payload'); base=[]; time=[];
                if validTrial, base=ctx.y; time=ctx.t; end
                taskShared{r}=struct('motion',file.payload,'baseline',base,'time',time,'trial',trialName);
            end
            for m=1:15
                index=index+1; prepared=cell(1,numel(uniquePP)); ppSeconds=zeros(1,numel(uniquePP));
                detectors=cell(1,numel(uniqueDetector));
                if validTrial
                    c=tbench.study.channel(ctx,m);
                    guard=tbench.study.guardView(c,cfg);
                else
                    c=struct('t',0,'fs',1,'envelope',NaN,'mvc',subject.mvc(m),'valid',false);
                    guard=struct('guardPeakHz',NaN,'guardStatistic',NaN,'guardAccepted',false,'guardOutsidePrimary',false);
                end
                for k=find(missing)
                    meta=struct('subject',subjects(s).id,'task',task,'trial',trialName,'muscle',m, ...
                        'muscleName',cfg.muscleNames{m},'method',candidates(k).analysisId, ...
                        'pipelineId',candidates(k).id,'candidateId',candidates(k).candidateId, ...
                        'candidateHash',candidates(k).hash,'preprocessId',candidates(k).preprocessId, ...
                        'inputKind',candidates(k).inputKind,'rawSamplingHz',rawFs,'analysisSamplingHz',c.fs, ...
                        'detectorRawFs',NaN,'runtimeSeconds',NaN,'preprocessSeconds',NaN,'inputValidFraction',0);
                    meta.guardPeakHz=guard.guardPeakHz; meta.guardStatistic=guard.guardStatistic;
                    meta.guardAccepted=guard.guardAccepted; meta.guardOutsidePrimary=guard.guardOutsidePrimary;
                    if validTrial
                        q=ppIndex(k); start=tic;
                        try
                            if isempty(prepared{q})
                                begin=tic; d=detectorIndex(k); p=candidates(k).preprocessParams;
                                if isempty(detectors{d}), detectors{d}=tbench.preprocess.detector(c,algorithmCfg,p); end
                                input=c; input.cachedDetector=detectors{d};
                                input.cachedDetectorParameters=[p.rawHighpassHz p.envelopePassHz p.envelopeStopHz];
                                prepared{q}=candidates(k).preprocess(input,algorithmCfg,candidates(k)); ppSeconds(q)=toc(begin);
                            end
                            meta.preprocessSeconds=ppSeconds(q); meta.detectorRawFs=prepared{q}.rawFs;
                            meta.inputValidFraction=mean(c.valid);
                            out=candidates(k).analyze(prepared{q},algorithmCfg,candidates(k));
                        catch err
                            if cfg.failFast, rethrow(err); end
                            out=tbench.emptyOutput(c); out.status='method_error'; out.message=err.message;
                        end
                        meta.runtimeSeconds=toc(start);
                    else
                        out=tbench.emptyOutput(c); out.status='invalid'; out.message=message;
                    end
                    [row,ev]=tbench.study.recordFeatures(out,c,cfg,meta);
                    if isempty(rows{k}), rows{k}=repmat(row,perSubject,1); end
                    rows{k}(index)=row; events{k}=append(events{k},ev);
                    if motionEnabled&&validTrial, taskPackets{k,r}=tbench.motion.packet(taskPackets{k,r},out,m,c); end
                end
            end
            if motionEnabled&&validTrial
                for k=find(missing), taskPackets{k,r}.records=rows{k}(index-14:index); end
            end
        end
        if motionEnabled
            for k=find(missing)
                detail=tbench.motion.task(taskPackets(k,:),taskShared,subjects(s),task,candidates(k),cfg);
                for name={'agreement','waveformPrediction','featureRows'}
                    key=name{1}; motionResults{k}.(key)=append(motionResults{k}.(key),detail.(key));
                end
            end
        end
        clear taskPackets taskShared;
    end
    inventories{s}=inv;
    if ~redLoaded, data.reduction=append(data.reduction,reduction); end
    for k=find(missing)
        payload=struct('key',keys{k},'sourceHash',sourceHashes{s},'codeHash',codeHash, ...
            'records',rows{k},'events',events{k},'inventory',inv,'reduction',reduction, ...
            'motionKey',motionKey,'motion',motionResults{k});
        if ~cfg.execution.cache||~isempty(sourceError), paths{k}=fullfile(cfg.workDir,'features',[keys{k} '.mat']); end
        tbench.study.atomicSave(paths{k},payload); data.recordFiles{s,k}=paths{k};
        data=addViews(data,rows{k},subjects(s),cfg,levels,s,k,N);
        if motionEnabled, data.motionViews{k}.rows=append(data.motionViews{k}.rows,motionResults{k}.featureRows); end
    end
end
level=find(strcmp(levels,cfg.clinical.featureLevel),1);
if isempty(level), error('tbench:FeatureLevel','Unsupported primary clinical level.'); end
for k=1:K, data.views{k}=data.levelViews{level,k}; end
data.inventory=concatenate(inventories); data.sourceHashes=sourceHashes;
end
function data=addViews(data,records,subject,cfg,levels,s,k,N)
for l=1:numel(levels)
    local=cfg; local.clinical.featureLevel=levels{l};
    [X,available,summary,names]=tbench.study.subjectFeatures(records,subject,local,levels{l});
    v=data.levelViews{l,k};
    if isempty(v.X), v.X=nan(N,size(X,2)); v.summary=repmat(summary,1,N); v.featureNames=names; end
    v.X(s,:)=X; v.available(s)=available; v.summary(s)=summary; data.levelViews{l,k}=v;
end
end
function a=append(a,b)
if isempty(b), return; end
if isempty(a), a=b(:); else, a=[a(:);b(:)]; end
end
function a=concatenate(C)
C=C(~cellfun('isempty',C)); a=struct([]);
if ~isempty(C), C=cellfun(@(v)v(:),C,'UniformOutput',false); a=vertcat(C{:}); end
end
