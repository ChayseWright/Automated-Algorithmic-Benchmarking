function result = task(packets,shared,subject,task,pipeline,cfg)
%TASK Compact M1-M4 rows from temporary packets, then release waveforms.
result=struct('agreement',struct([]),'waveformPrediction',struct([]),'featureRows',struct([]));
m=cfg.motion; m.sourceSynchronized=cfg.synchronized;
for r=1:numel(shared)
    if isempty(shared{r}), continue; end
    ref=shared{r}.motion; packet=packets{r};
    if isempty(ref.targets), continue; end
    for k=1:numel(ref.targets)
        target=ref.targets(k);
        for variant={'consensus','multitaper_only','oscillator_only'}
            witness=tbench.motion.referenceVariant(target,variant{1},ref,m);
            a=tbench.motion.scoreTrial(packet,ref,witness,m,pipeline.expectedDynamic); a.referenceVariant=variant{1};
            a=tag(a,subject,task,cfg.trials{r},pipeline); result.agreement=append(result.agreement,a);
        end
        features=emptyFeatures();
        if ~isempty(packet)&&~isempty(packet.records), features=tbench.motion.trialFeatures(packet.records,target.muscles); end
        row=struct('subject',subject.id,'task',task,'trial',cfg.trials{r},'target',target.name, ...
            'kind',target.kind,'unit',target.unit,'frame',target.frame,'referenceHash',ref.hash, ...
            'referenceValidFraction',target.validFraction,'inputCoverage',features.coverage, ...
            'meanPower',target.meanBandPower,'upperPower',target.upperWindowPower,'occupancy',target.occupancy, ...
            'F0',features.F0,'F1',features.F1,'F2',features.F2,'F3',features.F3,'R',features.R,'B0',features.B0);
        result.featureRows=append(result.featureRows,row);
        for lambda=cfg.motion.readoutRidges(:)'
            q=tbench.motion.waveformPrediction(packets,shared,r,target.name,m,lambda);
            q=tag(q,subject,task,cfg.trials{r},pipeline); result.waveformPrediction=append(result.waveformPrediction,q);
        end
    end
end
end
function x=emptyFeatures()
x=struct('F0',nan(1,2),'F1',nan(1,3),'F2',nan(1,6),'F3',nan(1,9),'R',nan(1,15),'B0',nan(1,2),'coverage',0);
end
function r=tag(r,s,task,trial,p)
r.subject=s.id; r.task=task; r.trial=trial; r.pipelineId=p.id; r.candidateId=p.candidateId; r.candidateHash=p.hash;
end
function a=append(a,b)
if isempty(b), return; end
if isempty(a), a=b(:); else, a=[a(:);b(:)]; end
end
