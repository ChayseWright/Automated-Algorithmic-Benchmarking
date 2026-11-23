function report = run_selftests()
%RUN_SELFTESTS Base MATLAB tests: no unittest or signal toolbox required.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root);
tests={@testTimeAndLayout,@testQuantile,@testSpectrum,@testDPSS,@testInputAliases, ...
    @testPhysicalScaling,@testCapabilities,@testReadout,@testFoldIsolation,@testIntegration};
report=struct('test',{},'passed',{},'message',{});
for k=1:numel(tests)
    name=func2str(tests{k});
    try
        tests{k}(); report(k)=struct('test',name,'passed',true,'message','');
        fprintf('PASS %s\n',name);
    catch err
        report(k)=struct('test',name,'passed',false,'message',err.message);
        fprintf('FAIL %s: %s\n',name,err.message);
    end
end
assert(all([report.passed]),'tbench:Selftests','One or more self-tests failed.');
end
function testTimeAndLayout()
cfg=tbench.defaultConfig(); fs=1000; t=(0:9999)'/fs; A=[zeros(numel(t),15) t];
tr=tbench.validateTrial(A,cfg); assert(abs(tr.fs-fs)<1e-8); assert(size(tr.x,2)==15);
B=A; B(:,16)=1e3*B(:,16); cfg.timeScale=.001; tr=tbench.validateTrial(B,cfg); assert(abs(tr.fs-fs)<1e-8);
cfg.timeScale=1; B=A; B(500,:)=[]; tr=tbench.validateTrial(B,cfg); assert(tr.missingSamples==1); assert(all(isnan(tr.x(500,:))));
bad=A; bad(400,16)=bad(399,16); assertThrows(@()tbench.validateTrial(bad,cfg));
assertThrows(@()tbench.validateTrial(A(:,1:15),cfg));
assert(~any(strcmp(cfg.tasks,'K5'))&&numel(cfg.tasks)==13);
end
function testQuantile()
q=tbench.quantile((1:5)',[.25 .75]); assert(isequal(size(q),[1 2]));
assert(isequal(q,[2 4])); assert(isequal(tbench.quantile(1:5,[.25 .75]),q));
assert(all(isnan(tbench.quantile([], [.25 .75]))));
end
function testSpectrum()
fs=100; t=(0:399)'/fs; y=cos(2*pi*6*t);
[P,f]=tbench.spectrum(y,fs,4,0,'hann'); assert(abs(trapz(f,P)-.5)<1e-6);
[~,i]=max(P); assert(abs(f(i)-6)<.2);
end
function testDPSS()
W=tbench.tapers(100,'dpss',1.5,2); assert(norm(W'*W-eye(2),'fro')<1e-10);
assert(isequal(W,tbench.tapers(100,'dpss',1.5,2)));
end
function testInputAliases()
data=struct('P4',struct('trial2',[zeros(100,15) (0:99)'/1000]));
a=tbench.loadSubject(struct('id','S01','EMGdata',data,'mvc',ones(1,15)));
assert(isfield(a.data.P4,'trial2'));
assertThrows(@()tbench.loadSubject(struct('id','S01','EMGData',data,'mvc',ones(1,14))));
file=[tempname '.mat']; cleanup=onCleanup(@()delete(file)); %#ok<NASGU>
EMGdata=data; MVC=zeros(100,15); save(file,'EMGdata','MVC');
a=tbench.loadSubject(struct('id','S01','file',file,'mvc',ones(1,15)));
assert(isequal(a.mvc,ones(1,15))); % An explicit reference overrides raw MVC traces.
end
function [ctx,cfg,tr]=toneContext()
cfg=tbench.defaultConfig(); fs=1000; t=(0:11999)'/fs;
x=(1+.65*cos(2*pi*6*t)).*cos(2*pi*90*t);
tr=tbench.validateTrial([repmat(x,1,15) t],cfg); ctx=tbench.prepare(tr,ones(1,15),cfg);
end
function testPhysicalScaling()
[ctx,cfg,tr]=toneContext(); c=tbench.channel(ctx,1); o=tbench.methodWelch(c,cfg,struct());
assert(o.detected&&abs(o.frequency-6)<.5);
cc=c; cc.envelope=2*c.envelope; cc.energy=4*c.energy; cc.y=2*c.y; cc.mvc=2;
p=tbench.methodWelch(cc,cfg,struct()); assert(abs(p.powerMean/o.powerMean-4)<1e-6);
tr.x=-tr.x; flipped=tbench.prepare(tr,ones(1,15),cfg);
v=ctx.valid&flipped.valid; difference=ctx.envelope(v)-flipped.envelope(v); assert(max(abs(difference))<1e-9);
assert(~isfield(c,'allMuscles')&&size(c.envelope,2)==1);
assert(size(c.raw,2)==1&&numel(c.rawTime)==size(tr.x,1)&&c.rawFs==tr.fs);
end
function testCapabilities()
[ctx,cfg]=toneContext(); c=tbench.channel(ctx,1);
o=tbench.validateOutput(tbench.methodWelch(c,cfg,struct()),c);
[r,ev]=tbench.summarize(o,c,cfg); assert(isnan(r.occupancy)&&isempty(ev));
o=tbench.validateOutput(tbench.methodLocal(c,cfg,struct('taper','hann')),c);
[r,~]=tbench.summarize(o,c,cfg); assert(isfinite(r.occupancy));
if r.occupancy>0, assert(abs(r.powerMean-r.occupancy*r.powerActive)<1e-10); end
broken=o; broken.time=[]; assertThrows(@()tbench.validateOutput(broken,c));
cc=c; cc.dynamicValid(:)=false; cc.y(:)=NaN;
insufficient=tbench.methodLocal(cc,cfg,struct('taper','hann'));
assert(strcmp(insufficient.status,'insufficient'));
tbench.validateOutput(insufficient,cc);
methods=tbench.defaultMethods();
for a=3:4
    candidate=methods(a).run(c,cfg,methods(a).params);
    candidate=tbench.validateOutput(candidate,c); [r,~]=tbench.summarize(candidate,c,cfg);
    assert(strcmp(candidate.status,'ok')&&isfinite(r.occupancy));
end
end
function testReadout()
X=ones(200,1); y=10*ones(200,1); fit=tbench.fitReadout(X,y,true(200,1),1);
p=tbench.predictReadout(fit,X); assert(max(abs(p-10))<1e-10);
end
function testFoldIsolation()
ids=arrayfun(@(s)sprintf('S%02d',s),1:25,'UniformOutput',false);
f=tbench.folds(ids,5,17); assert(numel(f)==25&&numel(unique({f.subject}))==25);
for k=1:5, assert(sum([f.fold]==k)==5); end
assert(isequal(f,tbench.folds(ids,5,17)));
end
function testIntegration()
folder=tempname; mkdir(folder); cleanup=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
cfg=tbench.defaultConfig(); cfg.expectedSubjects=1; cfg.tasks={'P1','K1'};
cfg.outputDir=folder; cfg.synchronized=true; cfg.verbose=false;
cfg.bench.calibration=false; cfg.bench.metrology=false; cfg.bench.network=false;
cfg.saveSignals=true; cfg.bootstrapSamples=20;
good=cfg.methods(1); bad=good; bad.name='intentional_failure'; bad.run=@intentionalFailure;
cfg.methods=[good bad]; subjects=tbench.demoSubjects(1,12,1000);
subjects.EMGData.P1=rmfield(subjects.EMGData.P1,'trial3');
s=tbench.run(subjects,cfg);
assert(numel(s.inventory)==6&&sum(strcmp({s.inventory.status},'invalid'))==1);
assert(numel(s.features)==180); % Both methods, 6 planned trials, 15 muscles.
assert(exist(fullfile(s.outputDir,'features.csv'),'file')==2);
assert(exist(fullfile(s.outputDir,'configuration.json'),'file')==2);
f=s.prediction(strcmp({s.prediction.method},'intentional_failure'));
assert(~isempty(f)&&all([f.predictorFailed])&&all([f.logScoreGain]==0));
missing=s.features(strcmp({s.features.status},'invalid')); assert(all(isnan([missing.powerMean])));
allFailed=s.repeatability(strcmp({s.repeatability.method},'intentional_failure'));
assert(~isempty(allFailed)&&all([allFailed.validTrials]==0)&&all([allFailed.missingTrials]==3));
cfg.calibration.nCases=1; cfg.calibration.duration=8; cfg.timeScale=.001; cfg.nominalFs=1111;
C=tbench.calibration(bad,cfg); assert(numel(C)==1&&strcmp(C.status,'method_error'));
tr=tbench.validateTrial(subjects.EMGData.P1.trial1,tbench.defaultConfig());
M=tbench.metrology(tr,ones(1,15),bad,tbench.defaultConfig());
assert(numel(M)==2&&all(strcmp({M.baseStatus},'method_error')));
end
function o=intentionalFailure(c,cfg,params) %#ok<INUSD,STOUT>
error('tbench:Intentional','Testing isolation of an algorithm failure.');
end
function assertThrows(f)
caught=false; try, f(); catch, caught=true; end
assert(caught,'Expected this invalid input to throw.');
end
