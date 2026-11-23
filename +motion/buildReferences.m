function result = buildReferences(subjects,cfg,codeHash,folder)
%BUILDREFERENCES Freeze one independent motion reference per acquired trial.
N=numel(subjects); J=numel(cfg.tasks); R=numel(cfg.trials); files=cell(N,J,R);
result=struct('subjectIds',{{subjects.id}},'files',{files},'subjectHashes',{cell(N,1)}, ...
    'referenceHashes',{cell(N,J,R)},'audit',struct([]),'cacheHits',0,'cacheMisses',0,'availableTrials',0,'referenceViews',[]);
if ~cfg.motion.enabled, return; end
for s=1:N
    source=struct(); identity=struct('status','missing_motion','source','','hash',''); sourceError='';
    try, [source,identity]=tbench.motion.loadSubject(subjects(s),cfg.motion);
    catch err, sourceError=err.message; identity.status='invalid_motion'; if cfg.failFast, rethrow(err); end
    end
    result.subjectHashes{s}=tbench.study.hash(identity);
    try, emg=tbench.loadSubject(subjects(s));
    catch, emg=struct('data',struct());
    end
    for j=1:J
        for r=1:R
            task=cfg.tasks{j}; trialName=cfg.trials{r}; status='missing_motion'; message=sourceError;
            ref=struct('schemaVersion',1,'status',status,'message',message,'targets',struct([]), ...
                'clock',tbench.motion.clock(struct(),[0;1],0,cfg.motion),'policyHash',tbench.motion.policyHash(cfg.motion));
            key=tbench.study.hash(struct('source',identity,'task',task,'trial',trialName,'code',codeHash,'policy',cfg.motion));
            try
                if ~strcmp(identity.status,'ok')||~isfield(source,task)||~isfield(source.(task),trialName)
                    error('tbench:MissingMotion','No canonical motion trial for %s.%s.',task,trialName);
                end
                if ~isfield(emg.data,task)||~isfield(emg.data.(task),trialName)
                    error('tbench:MissingEMGClock','No corresponding EMG clock for %s.%s.',task,trialName);
                end
                clockTrial=tbench.validateTrial(emg.data.(task).(trialName),cfg);
                clock=struct('timeOrigin',clockTrial.timeOrigin,'duration',clockTrial.duration,'samplingHz',clockTrial.fs);
                key=tbench.study.hash(struct('motion',key,'emgClock',clock));
                cache=fullfile(cfg.cacheDir,'motion_reference',[key '.mat']); hit=false;
                if cfg.execution.cache&&exist(cache,'file')==2
                    try
                        loaded=load(cache,'payload');
                        if strcmp(loaded.payload.key,key), ref=loaded.payload.reference; hit=true; end
                    catch
                    end
                end
                if hit, result.cacheHits=result.cacheHits+1;
                else
                    ref=tbench.motion.reference(source.(task).(trialName),clock,cfg.motion,task);
                    result.cacheMisses=result.cacheMisses+1;
                    if cfg.execution.cache, tbench.study.atomicSave(cache,struct('key',key,'reference',ref)); end
                end
                status='ok'; result.availableTrials=result.availableTrials+1;
            catch err
                status=err.identifier; message=err.message; ref.status=status; ref.message=message;
                if cfg.failFast&&~any(strcmp(err.identifier,{'tbench:MissingMotion','tbench:MissingEMGClock'})), rethrow(err); end
            end
            ref.hash=key; ref.subject=subjects(s).id; ref.task=task; ref.trial=trialName; ref.source=identity;
            safe=regexprep(subjects(s).id,'[^A-Za-z0-9_-]','_');
            path=fullfile(folder,'motion','references',safe,[task '_' trialName '.mat']);
            tbench.study.atomicSave(path,ref); result.files{s,j,r}=path;
            result.referenceHashes{s,j,r}=key;
            if isempty(ref.targets)
                row=auditRow(ref,subjects(s),task,trialName,status,message,[]);
                result.audit=append(result.audit,row);
            else
                for k=1:numel(ref.targets)
                    row=auditRow(ref,subjects(s),task,trialName,status,message,ref.targets(k));
                    result.audit=append(result.audit,row);
                end
            end
        end
    end
    result.subjectHashes{s}=tbench.study.hash(struct('source',identity,'references',{result.referenceHashes(s,:,:)}));
end
if cfg.motion.required&&result.availableTrials==0
    error('tbench:MotionRequired','No matched motion trials were imported. See the canonical MotionData schema and exact subject-ID join.');
end
result.referenceViews=tbench.motion.referenceView(result,subjects,cfg);
end
function r=auditRow(ref,subject,task,trial,status,message,target)
r=struct('subject',subject.id,'group',subject.group,'phenotype',subject.phenotype,'task',task,'trial',trial, ...
    'target','','kind','','unit','','status',status,'message',message,'referenceHash',ref.hash, ...
    'sourceHash',ref.source.hash,'aligned',ref.clock.aligned,'phaseReliable',ref.clock.phaseReliable, ...
    'timingReliable',ref.clock.timingReliable,'timingUncertaintySeconds',ref.clock.uncertaintySeconds, ...
    'clockOffsetSeconds',ref.clock.offsetSeconds,'clockRateRatio',ref.clock.rateRatio, ...
    'validFraction',0,'captureFraction',0,'recordedSeconds',0,'validSeconds',0,'nominalSeconds',NaN, ...
    'stateCoverage',0,'positiveSeconds',0,'negativeSeconds',0,'indeterminateSeconds',NaN, ...
    'meanBandPower',NaN,'limitRMS',NaN,'negativeReferenceCalibrated',false,'kineticReviewed',false);
if ~isempty(target)
    r.target=target.name; r.kind=target.kind; r.unit=target.unit;
    for key={'validFraction','captureFraction','recordedSeconds','validSeconds','nominalSeconds', ...
            'stateCoverage','positiveSeconds','negativeSeconds','indeterminateSeconds', ...
            'meanBandPower','limitRMS','negativeReferenceCalibrated','kineticReviewed'}
        r.(key{1})=target.(key{1});
    end
end
end
function a=append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
