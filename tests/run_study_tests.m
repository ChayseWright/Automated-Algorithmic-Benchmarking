function report = run_study_tests(includeEndToEnd,includeFullCatalog)
%RUN_STUDY_TESTS Native software-3 tests. Execute in MATLAB, not in Python.
if nargin<1, includeEndToEnd=false; end
if nargin<2, includeFullCatalog=false; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
tests={@testTransfer,@testRawFM,@testClockAndGaps,@testCatalogAndPairing, ...
    @testThresholdAndMetrics,@testTrainingOnlyClinical,@testSelection, ...
    @testCacheIdentity,@testCSV,@testFrozenCompression,@testRunNumbers,@testInstallation,@testAdapters,@testRepeatability, ...
    @testParameterSearch,@testCrossPipelines,@testFullRecordSubspace,@testPartialMerge,@testIRASAWarp};
if includeEndToEnd, tests{end+1}=@testEndToEnd; tests{end+1}=@testAdaptiveStudy; end
if includeFullCatalog, tests{end+1}=@testFullCatalog; end
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
assert(all([report.passed]),'tbench:StudyTests','One or more native study tests failed.');
end
function testTransfer()
cfg=tbench.studyConfig();
for fs=[1000 1111 2500]
    [h,a]=tbench.preprocess.designAA(fs,30,45,90,.05);
    assert(a.passErrorDB<=.05&&a.stopAttenuationDB>=90);
    assert(max(abs(h-flipud(h)))<1e-13);
    t=(0:12*fs-1)'/fs;
    for f=[3 6 9 18 27]
        [y,u]=tbench.preprocess.rateChange(cos(2*pi*f*t),fs,100,30,45,cfg.reduction);
        good=isfinite(y); X=[ones(size(u)) cos(2*pi*f*u) sin(2*pi*f*u)];
        b=pinv(X(good,:))*y(good); gain=norm(b(2:3)); phase=atan2(-b(3),b(2));
        assert(abs(gain-1)<.006&&abs(phase)<1e-3);
    end
    [y,u]=tbench.preprocess.rateChange(cos(2*pi*90*t),fs,100,30,45,cfg.reduction);
    good=isfinite(y); X=[ones(size(u)) cos(2*pi*90*u) sin(2*pi*90*u)];
    b=pinv(X(good,:))*y(good); assert(norm(b(2:3))<5e-5);
end
[y,u]=tbench.preprocess.rateChange(cos(2*pi*90*(0:11999)'/1000),1000,1000,30,45,cfg.reduction);
good=isfinite(y); X=[ones(size(u)) cos(2*pi*90*u) sin(2*pi*90*u)];
b=pinv(X(good,:))*y(good); assert(norm(b(2:3))<5e-5); % Full-rate comparison still applies the same observer.
end
function testRawFM()
cfg=tbench.studyConfig(); fs=2500; t=(0:12*fs-1)'/fs;
x=cos(2*pi*150*t+(10/6)*sin(2*pi*6*t));
[y,u]=tbench.preprocess.rateChange(x,fs,1250,450,537.5,cfg.reduction);
good=isfinite(y); expected=cos(2*pi*150*u+(10/6)*sin(2*pi*6*u));
assert(sqrt(mean((y(good)-expected(good)).^2))<1e-3);
end
function testClockAndGaps()
cfg=tbench.studyConfig(); fs=1111; t=(0:12*fs-1)'/fs;
x=(1+.6*cos(2*pi*6*t)).*cos(2*pi*110*t);
A=[repmat(x,1,15) t+42]; A(5000,:)=[];
trial=tbench.validateTrial(A,cfg); assert(trial.missingSamples==1&&all(isnan(trial.x(5000,:))));
ctx=tbench.preprocess.reference(trial,ones(1,15),cfg);
assert(abs(ctx.timeOrigin-42)<1e-10&&max(abs(diff(ctx.t)-.01))<1e-10);
gapTime=(5000-1)/trial.fs; assert(~any(ctx.valid(abs(ctx.t-gapTime)<.10,1)));
c=tbench.study.channel(ctx,1); assert(size(c.raw,2)==1&&~isfield(c,'group')&&~isfield(c,'severity_score'));
rails=cfg; rails.knownRails=[-1 1]; B=[zeros(numel(t),15) t]; B(100,1)=1;
masked=tbench.validateTrial(B,rails); assert(isnan(masked.x(100,1))&&isfinite(masked.x(99,1)));
p=tbench.pipelineRegistry('hann_tracker'); out=p.analyze(p.preprocess(c,cfg,p),cfg,p);
assert(~any(out.active(~out.valid)));
end
function testCatalogAndPairing()
cfg=tbench.studyConfig(); C=tbench.methodCatalog();
assert(numel(C)==112&&sum(strcmp({C.readiness},'implemented'))==58);
P=tbench.pipelineRegistry('implemented');
for k=1:numel(P), tbench.study.validatePair(P(k)); end
p=tbench.pipelineRegistry('cyclic_coherence'); p.inputKind='envelope_time';
mustThrow(@()tbench.study.validatePair(p),'tbench:Pair');
unavailable=find(~strcmp({C.readiness},'implemented'),1);
mustThrow(@()tbench.methodProposed(struct(),cfg,struct(),C(unavailable).id),'tbench:UnavailableModule');
K=tbench.study.candidates(cfg.pipelines,cfg); L=tbench.study.candidates(cfg.pipelines,cfg); assert(isequal({K.hash},{L.hash}));
assert(all(arrayfun(@(p)p.preprocessParams.envelopePassHz>=27,K)));
end
function testThresholdAndMetrics()
[tau,audit]=tbench.benchmark.controlThreshold([.1 .2 .2 .4 .9],true(1,5),repmat({'control'},1,5),.2);
assert(tau==.4&&audit.allowedErrors==1&&audit.empiricalFPR==.2);
[tau,audit]=tbench.benchmark.controlThreshold(.2,false,{'control'},.1);
assert(isinf(tau)&&audit.controlsUsed==0&&~audit.populationGuarantee);
cfg=tbench.studyConfig(); S=clinicalSubjects(6,6); X=[(1:12)' mod((1:12)',3) ones(12,1)];
model=tbench.benchmark.fitClinical(X,S,true(12,1),cfg,10);
P=tbench.benchmark.predictClinical(model,X,S,true(12,1),cfg,.5);
for k=1:numel(P)
    et=strcmp(P(k).group,'ET'); P(k).score=.1+.8*et; P(k).status='ok';
    P(k).severityPrediction=P(k).trueSeverityScore;
    g=P(k).trueSeverityGroup; P(k).ordinalP1=double(g==1); P(k).ordinalP2=double(g==2); P(k).ordinalP3=double(g==3);
end
P=tbench.benchmark.applyThreshold(P,.5,cfg); m=tbench.benchmark.metrics(P);
assert(m.controlFPR==0&&m.mildCapture==1&&m.severeCapture==1&&m.AUCavailable==1);
assert(abs(m.severitySpearman-1)<1e-12&&abs(m.severityKendallTauB-1)<1e-12&&m.severityRPS==0);
P(7).status='indeterminate'; P(7).decision=NaN; m=tbench.benchmark.metrics(P);
assert(m.mildCaptured<m.mildSubjects&&m.mildIndeterminate==1);
end
function testTrainingOnlyClinical()
cfg=tbench.studyConfig(); S=clinicalSubjects(6,6);
X=[(.1:.1:1.2)' (1:12)' zeros(12,1)]; train=[true(1,10) false false];
model=tbench.benchmark.fitClinical(X(train,:),S(train),true(10,1),cfg,10);
assert(all(ismember(model.trainingSubjects,{S(train).id})));
assert(model.featureTransform.scales(1)==median(X(1:6,1)));
before=tbench.benchmark.predictClinical(model,X(~train,:),S(~train),true(2,1),cfg,.5);
T=S; T(11).group='control'; T(11).severity_score=1e9; T(12).severity_score=-1e9;
after=tbench.benchmark.predictClinical(model,X(~train,:),T(~train),true(2,1),cfg,.5);
assert(isequal([before.score],[after.score]));
assert(isequal(model.severityCuts,tbench.quantile([S(7:10).severity_score],[1/3 2/3])));
plan=tbench.study.stratifiedFolds(clinicalSubjects(10,25),5,11);
assert(numel(unique({plan.subject}))==35);
for k=1:5
    assert(sum([plan.fold]==k&strcmp({plan.group},'control'))==2);
    assert(sum([plan.fold]==k&strcmp({plan.group},'ET'))==5);
end
end
function testSelection()
cfg=tbench.studyConfig(); m=struct('mildCapture',.9,'severeCapture',.5,'balancedLogLoss',.3);
s=repmat(struct('eligible',true,'violation',0,'metrics',m),1,3);
s(1).eligible=false; s(1).violation=.3;
s(2).metrics.mildCapture=.7; s(2).metrics.severeCapture=1;
s(3).metrics.mildCapture=.8; s(3).metrics.severeCapture=1;
[k,status]=tbench.study.chooseCandidate(s,cfg); assert(k==3&&strcmp(status,'constraints_met'));
[s.eligible]=deal(false); s(2).violation=.1; s(3).violation=.2;
[k,status]=tbench.study.chooseCandidate(s,cfg); assert(k==2&&strcmp(status,'constraints_unmet'));
end
function testCacheIdentity()
cfg=tbench.studyConfig(); cfg.tasks={'P1','K1'}; cfg.trials={'trial1'}; cfg.expectedDuration=5;
assert(~strcmp(tbench.study.hash(NaN),tbench.study.hash(Inf)));
assert(~strcmp(tbench.study.hash([1 2]),tbench.study.hash([1;2])));
S=tbench.demoStudySubjects(cfg,1,1); original=tbench.study.sourceHash(S(1));
T=S(1); T.group='ET'; T.severity_score=999; assert(strcmp(original,tbench.study.sourceHash(T)));
T.mvc(1)=3; assert(~strcmp(original,tbench.study.sourceHash(T)));
T=S(1); T.EMGData.P1.trial1(100,1)=T.EMGData.P1.trial1(100,1)+1;
assert(~strcmp(original,tbench.study.sourceHash(T)));
other=cfg; other.clinical.controlFPR=.2; other.outputDir=tempdir;
assert(strcmp(tbench.study.hash(tbench.study.processingConfig(cfg)),tbench.study.hash(tbench.study.processingConfig(other))));
other.analysisFs=200;
assert(~strcmp(tbench.study.hash(tbench.study.processingConfig(cfg)),tbench.study.hash(tbench.study.processingConfig(other))));
end
function testCSV()
file=[tempname '.csv']; cleanup=onCleanup(@()delete(file)); %#ok<NASGU>
fid=fopen(file,'w'); fprintf(fid,'subject,group,severity_score,phenotype\r\n"ET,01",ET,12,"ET-plus"\r\nC01,control,,unknown\r\n'); fclose(fid);
[h,r]=tbench.study.readCSV(file); assert(numel(h)==4&&strcmp(r{1}{1},'ET,01')&&isempty(r{2}{3}));
S=struct('id',{'C01','ET,01'}); S=tbench.attachClinical(S,file);
assert(isnan(S(1).severity_score)&&S(2).severity_score==12&&strcmp(S(2).phenotype,'ET-plus'));
end
function testFrozenCompression()
cfg=tbench.studyConfig(); S=clinicalSubjects(6,6); train=[true(1,10) false false]; test=~train;
rich=struct('X',repmat((1:12)',1,15)/10,'available',true(12,1), ...
    'featureNames',{{'meanPower','upperRecordPower','detectionYield','occupancy','activeWindowPower','hasActivity', ...
    'episodePower','episodeDuration','episodeRate','localPower50','localPower90','harmonicRatio','frequencySpread','frequencyChange','variance'}});
compressed=struct('X',rich.X(:,1:3),'available',rich.available,'featureNames',{rich.featureNames(1:3)});
local=cfg; local.clinical.featureLevel='R';
model=tbench.benchmark.fitClinical(rich.X(train,:),S(train),rich.available(train),local,10);
P=tbench.study.frozenCompression(rich,compressed,S,train,test,cfg,model,.5);
assert(numel(P)==2&&all(isfinite([P.score])));
end
function testRunNumbers()
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
[a,id1,n1]=tbench.study.newRun(folder,'test'); [b,id2,n2]=tbench.study.newRun(folder,'test');
assert(n2==n1+1&&~strcmp(a,b)&&~strcmp(id1,id2));
end
function testInstallation()
cfg=tbench.studyConfig(); R=tbench.validateInstallation(cfg); assert(numel(R)==6&&all(strcmp({R.status},'passed')));
end
function testAdapters()
cfg=tbench.studyConfig(); cfg.injection.muscles=1; cfg.injection.frequencies=6;
cfg.temporal.targets=1; fs=1000; t=(0:13999)'/fs;
raw=(1+.08*cos(2*pi*6*t)).*cos(2*pi*110*t);
trial=tbench.validateTrial([repmat(raw,1,15) t],cfg);
ctx=tbench.preprocess.reference(trial,ones(1,15),cfg);
p=tbench.pipelineRegistry('welch_cfar'); outputs=cell(1,15);
for m=1:15
    c=tbench.study.channel(ctx,m); outputs{m}=p.analyze(p.preprocess(c,cfg,p),cfg,p);
end
I=tbench.study.controlInjection(trial,ones(1,15),p,cfg,outputs,ctx);
assert(numel(I)==2&&all([I.truthPower]>0)&&all(~[I.waveformComparisonAvailable]));
assert(all(isnan([I.incrementalProjectedPower])));
T=tbench.study.temporalProbe(trial,ones(1,15),p,cfg,ctx);
assert(~isempty(T)&&all(isfinite([T.logScoreGain])));
end
function testRepeatability()
cfg=tbench.studyConfig(); cfg.tasks={'P1'};
row=struct('subject','','method','test','task','P1','muscle',1,'trial','', ...
    'status','ok','powerMean',NaN,'mvcPowerMean',NaN,'frequencyHz',6,'detected',true);
F=repmat(row,1,12); k=0;
for s=1:4
    for r=1:3
        k=k+1; F(k).subject=sprintf('S%d',s); F(k).trial=cfg.trials{r};
        F(k).powerMean=s+.05*r; F(k).mvcPowerMean=100*F(k).powerMean;
    end
end
[~,I,~]=tbench.repeatability(F,cfg); J=I([I.muscle]==1);
assert(J.ICC_A_1>.9&&isfinite(J.SEM_absolute)&&abs(J.MDC95-1.96*sqrt(2)*J.SEM_absolute)<1e-12);
end
function testFullCatalog()
cfg=tbench.studyConfig(); cfg.pipelines=tbench.pipelineRegistry('implemented');
% Detection failures are reported as failures, not installation exceptions.
R=tbench.validateInstallation(cfg); assert(numel(R)==58&&all([R.interfacePassed]));
disp(R(~[R.cleanToneRecovered]));
end
function testEndToEnd()
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
cfg=tbench.studyConfig(); cfg.expectedSubjects=12; cfg.expectedDuration=14;
cfg.tasks={'P1','K1'}; cfg.trials={'trial1','trial2'}; cfg.outerFolds=3;
cfg.tuning.innerFolds=2; cfg.tuning.maxCandidates=2; cfg.tuning.ridgeValues=[1 10];
cfg.tuning.requireParameterCoverage=false; cfg.tuning.refinementRounds=0;
cfg.tuning.bootstrapSamples=20; cfg.bench.calibration=false; cfg.tuning.calibrationGate=false;
cfg.execution.computeExplanatory=false; cfg.bench.metrology=false; cfg.verbose=false;
cfg.execution.candidateBatchSize=2; % Exercise multi-batch extraction/merge.
cfg.clinical.minSubjectCoverage=.45; cfg.outputDir=folder; cfg.cacheDir=fullfile(folder,'cache');
cfg.pipelines=tbench.pipelineRegistry({'welch_cfar','hann_tracker'});
S=tbench.demoStudySubjects(cfg,6,6); S(12).EMGData.P1=rmfield(S(12).EMGData.P1,'trial2');
study=tbench.runStudy(S,cfg);
P=study.clinicalPredictions(strcmp({study.clinicalPredictions.level},'F1'));
assert(numel(P)==24&&numel(study.inventory)==48&&sum(strcmp({study.inventory.status},'invalid'))==1);
for a=1:2
    Q=P(strcmp({P.pipelineId},cfg.pipelines(a).id)); assert(numel(unique({Q.subject}))==12);
end
for k=1:numel(study.tuningProvenance)
    q=study.tuningProvenance(k); train=strsplit(q.trainingSubjects); test=strsplit(q.validationSubjects);
    assert(isempty(intersect(train,test)));
    outerTest={study.outerPlan([study.outerPlan.fold]==q.outerFold).subject};
    assert(isempty(intersect([train test],outerTest)));
end
assert(exist(fullfile(study.outputDir,'selected_parameters.csv'),'file')==2);
assert(exist(fullfile(study.outputDir,'rhythmSummary.csv'),'file')==2);
files=dir(fullfile(study.outputDir,'selected','**','*.mat')); assert(~isempty(files));
p=tbench.loadSelected(fullfile(files(1).folder,files(1).name));
assert(numel(p.outputs{1}.time)==numel(p.time)&&isa(p.envelope,'double'));
cfg.execution.saveSelectedSignals=false; second=tbench.runStudy(S,cfg);
assert(second.runNumber==study.runNumber+1&&second.cacheHits==12*numel(study.candidateCatalog)&&second.cacheMisses==0);
assert(isequal({study.selection.candidateHash},{second.selection.candidateHash}));
end
function testParameterSearch()
cfg=tbench.studyConfig(); p=tbench.pipelineRegistry('dpss_tracker');
K=tbench.study.candidates(p,cfg); assert(numel(K)<=cfg.tuning.maxCandidates);
assert(all(arrayfun(@(q)q.analysisParams.nTapers<=2*q.analysisParams.NW-1,K)));
[pp,ap]=tbench.study.parameterSpace(p,cfg);
for module={'preprocessParams','analysisParams'}
    space=pp; if strcmp(module{1},'analysisParams'), space=ap; end
    names=fieldnames(space);
    for j=1:numel(names)
        key=names{j}; choices=space.(key); if ~iscell(choices), choices=num2cell(choices); end
        for v=1:numel(choices)
            assert(any(arrayfun(@(q)isequal(q.(module{1}).(key),choices{v}),K)));
        end
    end
end
cfg.tuning.spaces=struct('analysisId','welch_cfar', ...
    'preprocess',struct('rawHighpassHz',[20 35]),'analysis',struct('threshold',[2 4]));
cfg.tuning.searchMode='grid'; A=tbench.study.candidates(tbench.pipelineRegistry('welch_cfar'),cfg);
assert(numel(A)==4&&numel(unique({A.hash}))==4);
cfg.tuning.searchMode='coverage'; cfg.tuning.maxCandidates=2; cfg.tuning.requireParameterCoverage=false;
B=tbench.study.candidates(tbench.pipelineRegistry('welch_cfar'),cfg);
N=tbench.study.coordinateNeighbors(B(1),cfg,{B.hash},42); assert(~isempty(N));
assert(isempty(intersect({B.hash},{N.hash})));
other=tbench.study.coordinateNeighbors(B(1),cfg,{B.hash},42); assert(isequal({N.hash},{other.hash}));
end
function testCrossPipelines()
cfg=tbench.studyConfig(); P=tbench.pipelineRegistry({'hann_tracker','dpss_tracker'});
[X,rejected]=tbench.crossPipelines(P,P,cfg); assert(numel(X)==4&&isempty(rejected));
assert(numel(unique({X.id}))==4&&all(arrayfun(@(p)p.preprocessParams.freezeTransform,X)));
[~,ap]=tbench.study.parameterSpace(X(2),cfg); assert(~isfield(ap,'NW'));
P=tbench.pipelineRegistry({'cwt_morse','synchrosqueezing'});
if numel(P)==2
    [X,rejected]=tbench.crossPipelines(P,P,cfg); assert(numel(X)==2&&numel(rejected)==2);
end
end
function testFullRecordSubspace()
state=rng; cleanup=onCleanup(@()rng(state)); rng(77,'twister'); %#ok<NASGU>
cfg=tbench.studyConfig(); fs=100; t=(0:2999)'/fs;
y=.03*randn(size(t))+double(t>=20&t<28).*cos(2*pi*6*t); valid=true(size(t));
c=struct('t',t,'fs',fs,'y',y,'envelope',y+2,'valid',valid,'dynamicValid',valid);
for kind={'music','esprit','pencil'}
    out=tbench.methodSubspace(c,cfg,struct('embedding',64,'signalRank',2,'lineFThreshold',4),kind{1});
    assert(out.detected&&abs(out.frequency-6)<.10&&~out.audit.timeTruncated);
    assert(out.audit.inferenceSamples==numel(t)&&out.audit.inferenceCoverageFraction==1);
end
c.y(1001:1020)=NaN; c.valid(1001:1020)=false; c.dynamicValid=c.valid;
out=tbench.methodSubspace(c,cfg,struct('embedding',64,'signalRank',2),'music');
assert(out.audit.inferenceCoverageFraction==1&&out.audit.inferenceSamples==sum(c.valid));
% Tiny finite-iteration atomic contract probe, not a convergence claim.
c=struct('t',(0:159)'/fs,'fs',fs,'y',cos(2*pi*6*(0:159)'/fs), ...
    'valid',true(160,1),'dynamicValid',true(160,1)); c.envelope=c.y+2;
out=tbench.methodSubspace(c,cfg,struct('atomicSamples',64,'iterations',2),'atomic');
assert(out.audit.inferenceCoverageFraction==1&&out.audit.blocks(end).stop>=1.6-1e-12);
end
function testPartialMerge()
cfg=tbench.studyConfig(); data=struct('recordFiles',{cell(4,1)},'sourceHashes',{{'a';'b';'c';'d'}}, ...
    'levelViews',{cell(6,1)},'views',{cell(1,1)},'cacheHits',0,'cacheMisses',0);
addition=struct('recordFiles',{{'one';'two'}},'sourceHashes',{{'b';'d'}}, ...
    'levelViews',{cell(6,1)},'cacheHits',1,'cacheMisses',1);
for l=1:6
    addition.levelViews{l,1}=struct('X',[10 20;30 40],'available',[true;false], ...
        'summary',repmat(struct('subject','test'),1,2),'featureNames',{{'x','y'}});
end
merged=tbench.study.mergeCandidateData(data,addition,[2 4],2,cfg);
assert(isequal(merged.views{2}.X([2 4],:),[10 20;30 40]));
assert(all(all(isnan(merged.views{2}.X([1 3],:))))&&~any(merged.views{2}.available([1 3])));
addition.sourceHashes{1}='CHANGED';
mustThrow(@()tbench.study.mergeCandidateData(data,addition,[2 4],2,cfg),'tbench:InputChanged');
end
function testIRASAWarp()
cfg=tbench.studyConfig(); fs=100; t=(0:1999)'/fs;
for factor=[1/1.9 1.9]
    y=tbench.preprocess.irasaWarp(cos(2*pi*6*t),factor,fs,cfg.reduction);
    u=(0:numel(y)-1)'/fs; f=6/factor; good=isfinite(y);
    X=[ones(size(u)) cos(2*pi*f*u) sin(2*pi*f*u)]; b=pinv(X(good,:))*y(good);
    assert(abs(norm(b(2:3))-1)<.005);
end
end
function testAdaptiveStudy()
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
cfg=tbench.studyConfig(); cfg.expectedSubjects=12; cfg.expectedDuration=14;
cfg.tasks={'P1','K1'}; cfg.trials={'trial1'}; cfg.outerFolds=3; cfg.tuning.innerFolds=2;
cfg.tuning.maxCandidates=2; cfg.tuning.requireParameterCoverage=false;
cfg.tuning.refinementRounds=1; cfg.tuning.refinementCandidates=2; cfg.tuning.refinementSeeds=1;
cfg.tuning.ridgeValues=10; cfg.tuning.bootstrapSamples=10;
cfg.tuning.spaces=struct('analysisId','welch_cfar','preprocess',struct('rawHighpassHz',[20 35]), ...
    'analysis',struct('threshold',[2 4 6]));
cfg.bench.calibration=false; cfg.tuning.calibrationGate=false; cfg.bench.metrology=false;
cfg.execution.computeExplanatory=false; cfg.execution.computeCompression=false;
cfg.execution.saveSelectedSignals=false; cfg.execution.validateModules=false; cfg.verbose=false;
cfg.outputDir=folder; cfg.cacheDir=fullfile(folder,'cache'); cfg.pipelines=tbench.pipelineRegistry('welch_cfar');
S=tbench.demoStudySubjects(cfg,6,6); study=tbench.runStudy(S,cfg);
assert(numel(study.candidateCatalog)>2&&any([study.selectionScores.searchRound]>0));
for k=1:numel(study.optimizationTrace)
    row=study.optimizationTrace(k); ids=strsplit(row.trainingSubjects);
    if row.outerFold<=cfg.outerFolds
        test={study.outerPlan([study.outerPlan.fold]==row.outerFold).subject}; assert(isempty(intersect(ids,test)));
    else, assert(numel(ids)==12&&strcmp(row.scope,'all_training_deployment'));
    end
end
assert(numel(study.clinicalPredictions)==12);
assert(exist(fullfile(study.outputDir,'parameter_spaces.csv'),'file')==2);
end
function S=clinicalSubjects(nc,ne)
blank=struct('id','','group','control','severity_score',NaN,'severity_kind','synthetic_ordered', ...
    'severity_group','','phenotype','synthetic');
S=repmat(blank,1,nc+ne);
for k=1:numel(S)
    S(k).id=sprintf('TEST_%02d',k);
    if k>nc, S(k).group='ET'; S(k).severity_score=k-nc; end
end
end
function mustThrow(f,identifier)
caught=false;
try, f(); catch err, caught=strcmp(err.identifier,identifier); end
assert(caught,'Expected exception %s.',identifier);
end
