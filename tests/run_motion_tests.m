function report = run_motion_tests(includeEndToEnd)
%RUN_MOTION_TESTS Native MATLAB regression and synthetic integration tests.
% This suite must execute in MATLAB before interpreting patient results.
if nargin<1, includeEndToEnd=false; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
tests={@testMotionConfig,@testClock,@testPhysicalSchema,@testPoseMasks,@testPolicyIdentity, ...
    @testMotionReference,@testCalibrationGate,@testSampling,@testEpisodeAccounting, ...
    @testCoherenceCapacity,@testNoTargetLeakage,@testMechanicalReadout,@testWaveformReadout, ...
    @testMotionSelection,@testPartialMotionMerge,@testUnknownReference};
if includeEndToEnd, tests{end+1}=@testReferenceCache; tests{end+1}=@testMotionStudy; end
report=struct('test',{},'passed',{},'message',{},'engine',{});
for k=1:numel(tests)
    try
        tests{k}(); report(k)=struct('test',func2str(tests{k}),'passed',true,'message','','engine',version);
        fprintf('PASS %s\n',report(k).test);
    catch err
        report(k)=struct('test',func2str(tests{k}),'passed',false,'message',err.message,'engine',version);
        fprintf('FAIL %s: %s\n',report(k).test,err.message);
    end
end
assert(all([report.passed]),'tbench:MotionTests','One or more native motion tests failed.');
end
function testMotionConfig()
c=tbench.studyConfig(); c.motion.enabled=true; tbench.motion.validateConfig(c.motion,c);
m=c.motion; m.maxFs=20; mustThrow(@()tbench.motion.validateConfig(m,c),'tbench:MotionConfig');
m=c.motion; m.searchBand=[4 8]; mustThrow(@()tbench.motion.validateConfig(m,c),'tbench:MotionConfig');
m=c.motion; m.enabled=NaN; mustThrow(@()tbench.motion.validateConfig(m,c),'tbench:MotionConfig');
end
function testClock()
m=tbench.motionConfig(); time=(500:.01:530)';
s=struct('motionEvents',[500;530],'emgEvents',[42;42+30*1.0002],'uncertaintySeconds',.001);
c=tbench.motion.clock(s,time,42,m);
assert(c.aligned&&c.phaseReliable&&abs(c.rateRatio-1.0002)<1e-12);
assert(max(abs(c.offsetSeconds+c.rateRatio*s.motionEvents-s.emgEvents))<1e-10);
s.uncertaintySeconds=.006; c=tbench.motion.clock(s,time,42,m); assert(c.timingReliable&&~c.phaseReliable);
c=tbench.motion.clock(struct(),time,42,m); assert(~c.aligned&&~c.phaseReliable);
s=struct('motionEvents',500,'emgEvents',42,'uncertaintySeconds',.001);
c=tbench.motion.clock(s,time,42,m); assert(~c.aligned); s.clockPPMbound=100;
c=tbench.motion.clock(s,time,42,m); assert(c.aligned&&abs(c.uncertaintySeconds-.004)<1e-10);
s.clockPPMbound=-1; mustThrow(@()tbench.motion.clock(s,time,42,m),'tbench:MotionSync');
end
function testPhysicalSchema()
m=tbench.motionConfig(); r=motionRecord(100,12,6); r.targets.values=r.targets.values/1000; r.targets.unit='m';
a=tbench.motion.validate(r,struct('timeOrigin',0),m); assert(strcmp(a.targets.unit,'mm'));
assert(abs(a.targets.values(1,1)-1)<1e-10);
r.time(100)=[]; r.targets.values(100,:)=[];
a=tbench.motion.validate(r,struct('timeOrigin',0),m); assert(a.missingSamples==1&&all(isnan(a.targets.values(100,:))));
r=motionRecord(100,12,6); r.targets.valid=true(numel(r.time),1); r.targets.valid=double(r.targets.valid); r.targets.valid(5)=NaN;
mustThrow(@()tbench.motion.validate(r,struct('timeOrigin',0),m),'tbench:MotionQuality');
r=motionRecord(100,12,6); r.targets.kind='joint_angle'; r.targets.unit='deg';
r.targets.values=mod((1:numel(r.time))'+170,360)-180;
a=tbench.motion.validate(r,struct('timeOrigin',0),m); assert(max(abs(diff(a.targets.values)))<pi/180+1e-10);
end
function testPoseMasks()
t=(0:99)'/100; p=struct('name','thorax','position',zeros(100,3), ...
    'quaternion',repmat([1 0 0 0],100,1),'valid',true(100,1));
S=[p p]; S(2).name='hand'; S(2).position(:,1)=1; S(2).valid(20)=false;
r=tbench.motion.fromPoses(t,S,'thorax','hand','mm',{'thorax','hand'});
assert(all(r.targets(1).values(1,:)==[1 0 0])&&all(isnan(r.targets(1).values(20,:))));
assert(all(isnan(r.targets(2).values(20,:)))&&strcmp(r.targets(2).kind,'rotation'));
end
function testPolicyIdentity()
m=tbench.motionConfig(); a=tbench.motion.policyHash(m); n=m; n.readoutRidges=1;
assert(strcmp(a,tbench.motion.policyHash(n))); n=m; n.windowSeconds=4;
assert(~strcmp(a,tbench.motion.policyHash(n)));
end
function testMotionReference()
m=tbench.motionConfig(); r=motionRecord(333,16,6);
a=tbench.motion.reference(r,struct('timeOrigin',0),m,'P1'); q=a.targets(1);
assert(a.clock.phaseReliable&&q.fs==100&&abs(q.meanBandPower-.52)<.04);
assert(any(q.multitaperState==1)&&isfinite(median(q.frequencyTrack(q.multitaperState==1))));
assert(any(q.state==1));
assert(abs(median(q.frequencyTrack(q.multitaperState==1))-6)<.25);
assert(~any(q.state==-1)&&~q.negativeReferenceCalibrated);
assert(q.transfer.stopAttenuationDB>=90&&q.transfer.passErrorDB<=.05);
r.reviewed=false; b=tbench.motion.reference(r,struct('timeOrigin',0),m,'K1'); assert(~any(b.targets.state~=0));
r.targets.confounded=true(numel(r.time),1); b=tbench.motion.reference(r,struct('timeOrigin',0),m,'P1');
assert(all(b.targets.state==0));
r=motionRecord(100,12,6); b=tbench.motion.reference(r,struct('timeOrigin',0,'duration',30),m,'P1');
assert(b.targets.validFraction<.4&&b.targets.captureFraction<=.4+1e-12);
end
function testCalibrationGate()
m=tbench.motionConfig(); r=motionRecord(100,12,6);
r.targets.noiseSD=.01; r.targets.limitRMS=.1; r.targets.calibrationVerified=true;
r.targets.calibrationPolicyHash='stale_policy';
a=tbench.motion.reference(r,struct('timeOrigin',0),m,'P1'); assert(~a.targets.negativeReferenceCalibrated);
r.targets.calibrationPolicyHash=tbench.motion.policyHash(m);
a=tbench.motion.reference(r,struct('timeOrigin',0),m,'P1'); assert(a.targets.negativeReferenceCalibrated);
r.targets=rmfield(r.targets,'noiseSD'); a=tbench.motion.reference(r,struct('timeOrigin',0),m,'P1');
assert(~a.targets.negativeReferenceCalibrated);
end
function testSampling()
m=tbench.motionConfig(); fs=1000; t=(0:11999)'/fs;
for f=[3 6 9]
    target=struct('time',(0:1199)'/100+.0007,'fs',100);
    [y,good]=tbench.motion.align(t,cos(2*pi*f*t),target,m);
    assert(mean(good)>.8&&sqrt(mean((y(good)-cos(2*pi*f*target.time(good))).^2))<.006);
end
x=cos(2*pi*6*t); x(t>5&t<5.1)=NaN; [~,good]=tbench.motion.align(t,x,target,m);
assert(~any(good(abs(target.time-5.05)<.1)));
end
function testEpisodeAccounting()
[p,ref,target,m]=scoringFixture(); q=tbench.motion.scoreTrial(p,ref,target,m,true);
assert(abs(q.recall-.5)<1e-12&&abs(q.missedMotionBurden-.5)<1e-12);
assert(abs(q.frequencyWithinToleranceRecall-.5)<1e-12&&abs(q.trackingMissingFraction-.5)<1e-12);
assert(abs(q.FPseconds-1)<1e-12&&abs(q.unconfirmedPositiveSeconds-1)<1e-12);
assert(abs(q.falseAlarmSecondsPerMinute-15)<1e-12);
q=tbench.motion.scoreTrial(p,ref,target,m,false); assert(isnan(q.recall)&&~q.nativeTemporalAvailable);
q=tbench.motion.scoreTrial([],ref,target,m,true); assert(q.recall==0&&q.missedMotionBurden==1&&q.FNseconds==2);
target.muscles=[1 2]; p.frequency=[6 6]; p.frequencyTrack=repmat(p.frequencyTrack,1,2);
p.valid=repmat(p.valid,1,2); p.active=repmat(p.active,1,2);
q=tbench.motion.scoreTrial(p,ref,target,m,true); assert(isnan(q.recall)&&strcmp(q.status,'EMG_source_clock_unverified'));
m.sourceSynchronized=true; q=tbench.motion.scoreTrial(p,ref,target,m,true); assert(q.recall==.5);
ref.clock.timingReliable=false; q=tbench.motion.scoreTrial(p,ref,target,m,true); assert(isnan(q.recall));
end
function testCoherenceCapacity()
m=tbench.motionConfig(); m.coherence.minimumSegments=3; m.coherence.surrogates=0;
m.sourceSynchronized=true;
t=(0:2999)'/100; x=cos(2*pi*6*t); target=struct('name','test','unit','mm','fs',100, ...
    'time',t,'waveform',x,'valid',true(size(t)),'muscles',[1 2],'frequency',6);
p=struct('time',t,'waveform',[x x],'valid',true(numel(t),2),'hasWaveform',[true true]);
ref=struct('clock',struct('phaseReliable',true)); q=tbench.motion.coherence(p,ref,target,m);
assert(all([q.coherentFraction]>.99)&&all([q.inputRank]==1)&&all(isnan([q.excludedContribution])));
assert(all(isnan([q.maximumSearchP])));
target.muscles=1; q=tbench.motion.coherence(p,ref,target,m);
assert(q.linearModelConditionsMet&&q.isolatedContribution>.99&&q.measuredMotionPower>0);
end
function testNoTargetLeakage()
c=tbench.studyConfig(); c.motion.enabled=true; c.motionRuntime=struct('secretMotion',123);
a=tbench.motion.inferenceConfig(c);
assert(~any(isfield(a,{'motion','motionRuntime','motionProtocol','clinical','tuning'})));
[S,V,c]=readoutFixture(); train=[true(1,10) false false];
plan=tbench.study.stratifiedFolds(S(train),2,17);
before=tbench.motion.trainingScore(V,S,train,plan,c);
for k=1:numel(V.rows)
    if any(strcmp(V.rows(k).subject,{S(~train).id})), V.rows(k).meanPower=1e30; end
end
after=tbench.motion.trainingScore(V,S,train,plan,c);
assert(before.available&&isequaln(before,after));
end
function testMechanicalReadout()
[S,V,c]=readoutFixture(); training={S(1:10).id}; test={S(11:12).id};
[A,models]=tbench.motion.lossRows(V,training,test,c,'F1',1);
assert(numel(A)==4&&numel(models)==1&&all(isfinite([A.predictedValue])));
for k=1:numel(V.rows)
    if any(strcmp(V.rows(k).subject,test)), V.rows(k).meanPower=100*V.rows(k).meanPower; end
end
[B,other]=tbench.motion.lossRows(V,training,test,c,'F1',1);
assert(isequal([A.predictedValue],[B.predictedValue])&&strcmp(tbench.study.hash(models),tbench.study.hash(other)));
assert(any([A.transformedLoss]~=[B.transformedLoss]));
for k=1:numel(V.rows), if any(strcmp(V.rows(k).subject,test)), V.rows(k).meanPower=NaN; end; end
A=tbench.motion.lossRows(V,training,test,c,'F1',1); assert(isempty(A));
mustThrow(@()tbench.motion.lossRows(V,training,training,c,'F1',1),'tbench:MotionLeakage');
end
function testWaveformReadout()
m=tbench.motionConfig(); t=(0:1599)'/100; packets=cell(1,3); shared=packets;
for k=1:3
    x=cos(2*pi*6*t+.4*k); y=.7*cos(2*pi*6*t+.4*k-.3);
    target=struct('name','test','unit','mm','kind','translation','frame','thorax', ...
        'muscles',1,'fs',100,'time',t,'waveform',y,'valid',true(size(t)));
    packets{k}=struct('time',t,'waveform',x,'valid',true(size(t)),'hasWaveform',true,'records',struct([]));
    shared{k}=struct('time',t,'baseline',x,'trial',sprintf('trial%d',k), ...
        'motion',struct('targets',target,'clock',struct('phaseReliable',true)));
end
q=tbench.motion.waveformPrediction(packets,shared,3,'test',m,1);
assert(strcmp(q.status,'ok')&&q.gainOverNull>.95&&q.inputCoverage>.99);
assert(strcmp(q.trainingTrials,'trial1 trial2'));
shared{3}.motion.clock.phaseReliable=false;
q=tbench.motion.waveformPrediction(packets,shared,3,'test',m,1); assert(strcmp(q.status,'phase_alignment_unavailable'));
end
function testMotionSelection()
[S,V,c]=readoutFixture(); c.tuning.maxCandidates=1; c.tuning.requireParameterCoverage=false;
c.tuning.enabled=false; c.tuning.calibrationGate=false; c.bench.calibration=false;
p=tbench.study.candidates(tbench.pipelineRegistry('welch_cfar'),c); p=p(1);
view=struct('X',repmat((1:12)',1,3),'available',true(12,1)); training=[true(1,10) false false];
[winner,~]=tbench.study.tune(S,training,1,p,{view},struct([]),c,1,{V});
assert(winner.motion.available&&winner.motion.subjects>0); % Regression: ninth argument must not be discarded.
s=repmat(struct('eligible',true,'violation',0,'metrics',struct('mildCapture',.8,'severeCapture',1,'balancedLogLoss',.2), ...
    'motion',struct('gain',0)),1,2); s(2).motion.gain=.5;
[k,~]=tbench.study.chooseCandidate(s,c); assert(k==2);
s(2).metrics.mildCapture=.7; [k,~]=tbench.study.chooseCandidate(s,c); assert(k==1);
end
function testPartialMotionMerge()
c=tbench.studyConfig(); data=struct('recordFiles',{cell(2,1)},'sourceHashes',{{'a';'b'}}, ...
    'levelViews',{cell(6,1)},'views',{cell(1,1)},'cacheHits',0,'cacheMisses',0, ...
    'motionViews',{{struct('rows',struct('subject','S2','target','old'))}});
addition=struct('recordFiles',{{'two'}},'sourceHashes',{{'b'}},'levelViews',{cell(6,1)}, ...
    'motionViews',{{struct('rows',struct([]))}},'cacheHits',0,'cacheMisses',1);
for l=1:6
    addition.levelViews{l,1}=struct('X',[1 2],'available',true,'summary',struct('subject','S2'),'featureNames',{{'x','y'}});
end
merged=tbench.study.mergeCandidateData(data,addition,2,1,c); assert(isempty(merged.motionViews{1}.rows));
end
function testUnknownReference()
[p,ref,target,m]=scoringFixture(); target.state(:)=0; target.positiveSeconds=0; target.negativeSeconds=0;
target.events=struct([]); q=tbench.motion.scoreTrial(p,ref,target,m,true);
assert(q.FPseconds==0&&isnan(q.falseAlarmSecondsPerMinute)&&isnan(q.recall)&&q.unconfirmedPositiveSeconds>0);
end
function testReferenceCache()
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
c=tbench.studyConfig(); c.expectedDuration=12; c.tasks={'P1','K1'}; c.trials={'trial1'};
c.motion.enabled=true; c.motion.required=true; c.cacheDir=fullfile(folder,'cache');
S=tbench.demoMotionSubjects(c,1,0); A=tbench.motion.buildReferences(S,c,'test',fullfile(folder,'one'));
S.group='ET'; S.severity_score=999;
B=tbench.motion.buildReferences(S,c,'test',fullfile(folder,'two'));
assert(B.cacheHits==2&&isequal(A.subjectHashes,B.subjectHashes));
S.MotionData=struct(); mustThrow(@()tbench.motion.buildReferences(S,c,'test',fullfile(folder,'empty')),'tbench:MotionRequired');
end
function testMotionStudy()
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
c=tbench.studyConfig(); c.expectedSubjects=12; c.expectedDuration=16; c.tasks={'P1','K1'};
c.outerFolds=3; c.tuning.innerFolds=2; c.tuning.maxCandidates=1; c.tuning.requireParameterCoverage=false;
c.tuning.refinementRounds=0; c.tuning.ridgeValues=10; c.tuning.bootstrapSamples=10;
c.bench.calibration=false; c.tuning.calibrationGate=false; c.bench.metrology=false;
c.execution.computeExplanatory=false; c.execution.validateModules=false; c.verbose=false;
c.motion.enabled=true; c.motion.required=true; c.motion.calibration.enabled=false; c.motion.coherence.enabled=false;
c.synchronized=true;
c.motion.minimumFitRows=8; c.motion.minimumMotionSubjects=2; c.motion.readoutRidges=[1 10];
c.pipelines=tbench.pipelineRegistry({'welch_cfar','hann_tracker'}); c.outputDir=folder; c.cacheDir=fullfile(folder,'cache');
S=tbench.demoMotionSubjects(c,6,6); S(12).MotionData.P1=rmfield(S(12).MotionData.P1,'trial2');
study=tbench.runStudy(S,c); P=study.clinicalPredictions(strcmp({study.clinicalPredictions.level},'F1'));
assert(numel(P)==24&&numel(study.motion.clinicalPredictions)==12);
assert(numel(study.motion.referenceAudit)==72&&sum(strcmp({study.motion.referenceAudit.status},'tbench:MissingMotion'))==1);
assert(~isempty(study.motion.predictions)&&~isempty(study.motion.waveformPrediction));
for k=1:numel(study.motion.predictions)
    r=study.motion.predictions(k); assert(~any(strcmp(strsplit(r.trainingSubjects),r.subject)));
end
A=study.motion.agreement(strcmp({study.motion.agreement.pipelineId},c.pipelines(1).id));
assert(all(isnan([A.recall]))&&all(~[A.nativeTemporalAvailable]));
assert(exist(fullfile(study.outputDir,'motion_coverage.csv'),'file')==2);
assert(exist(fullfile(study.outputDir,'motion_agreementSummary.csv'),'file')==2);
files=dir(fullfile(study.outputDir,'selected','**','*.mat')); p=tbench.loadSelected(fullfile(files(1).folder,files(1).name));
assert(isfield(p,'motionReference')&&numel(p.outputs{1}.time)==numel(p.time));
end
function r=motionRecord(fs,T,f)
t=(0:round(fs*T)-1)'/fs;
target=struct('name','hand','kind','translation','unit','mm','frame','thorax', ...
    'values',[cos(2*pi*f*t)+.1*sin(2*pi*.3*t),.2*sin(2*pi*f*t),zeros(size(t))]);
r=struct('time',t,'targets',target,'reviewed',true, ...
    'sync',struct('verified',true,'offsetSeconds',0,'rateRatio',1,'uncertaintySeconds',0));
end
function [p,ref,target,m]=scoringFixture()
m=tbench.motionConfig(); t=(0:999)'/100; state=zeros(1000,1,'int8'); state(101:300)=1; state(601:1000)=-1;
target=struct('name','hand','kind','translation','unit','mm','time',t,'fs',100,'muscles',1, ...
    'state',state,'valid',true(size(t)),'frequency',6,'frequencyTrack',6*ones(size(t)), ...
    'positiveSeconds',2,'negativeSeconds',4,'indeterminateSeconds',4,'validFraction',1,'stateCoverage',.6, ...
    'negativeReferenceCalibrated',true,'powerSigned',ones(size(t)));
target.events=tbench.episodes(state==1,t,target.frequencyTrack,target.powerSigned,target.valid,0);
p=struct('time',t,'frequency',6,'frequencyTrack',6*ones(size(t)),'dynamic',true,'valid',true(size(t)),'active',false(size(t)));
p.active([101:200 401:500 601:700])=true;
ref=struct('clock',struct('uncertaintySeconds',0,'timingReliable',true));
end
function [S,V,c]=readoutFixture()
c=tbench.studyConfig(); c.tasks={'P1','K1'}; c.trials={'trial1'}; c.motion.enabled=true;
c.motion.minimumFitRows=2; c.motion.minimumMotionSubjects=1; c.motion.readoutRidges=[1 10]; c.tuning.innerFolds=2;
blank=struct('id','','group','control','severity_score',NaN,'severity_kind','synthetic','severity_group','','phenotype','synthetic');
S=repmat(blank,1,12); rows=struct([]);
for s=1:12
    S(s).id=sprintf('T%02d',s); if s>6, S(s).group='ET'; S(s).severity_score=s-6; end
    for j=1:2
        f=[s/10 s/5 .5]; row=struct('subject',S(s).id,'task',c.tasks{j},'trial','trial1', ...
            'target','hand','kind','translation','unit','mm','frame','thorax','referenceHash','synthetic', ...
            'referenceValidFraction',1,'inputCoverage',1,'meanPower',(s/10)^2, ...
            'B0',[1 1],'F0',f(1:2),'F1',f,'F2',[f zeros(1,3)],'F3',[f zeros(1,6)],'R',[f zeros(1,12)]);
        if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
    end
end
V=struct('rows',rows);
end
function mustThrow(f,identifier)
caught=false; try, f(); catch err, caught=strcmp(err.identifier,identifier); end
assert(caught,'Expected exception %s.',identifier);
end
