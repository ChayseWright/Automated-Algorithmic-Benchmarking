function r = waveformPrediction(packets,shared,held,targetName,m,lambda)
%WAVEFORMPREDICTION Two-repetition calibration, third repetition evaluation.
r=struct('target',targetName,'unit','','status','unavailable','ridge',lambda, ...
    'scope','person_calibrated_leave_one_repetition_out','trainingTrials','','baselineRidge',m.baselineRidge, ...
    'testRows',0,'fitRows',0,'inputCoverage',NaN,'nullMSE',NaN,'referenceMSE',NaN,'methodMSE',NaN, ...
    'gainOverNull',NaN,'gainOverReference',NaN,'logScoreGainOverReference',NaN, ...
    'observationVariance',NaN,'observationScaleScope','training_baseline_residual_isotropic_with_relative_floor', ...
    'frozenTemporalMeanMSE',NaN,'frozenTemporalMeanRetained',NaN, ...
    'compressionScope','trial_constant_F0_F1_F2_F3_Rsubject_and_frozen_temporal_mean', ...
    'F0waveformMSE',NaN,'F1waveformMSE',NaN,'F2waveformMSE',NaN,'F3waveformMSE',NaN,'RsubjectWaveformMSE',NaN, ...
    'F0waveformRetained',NaN,'F1waveformRetained',NaN,'F2waveformRetained',NaN,'F3waveformRetained',NaN,'RsubjectWaveformRetained',NaN);
if isempty(shared{held}), return; end
reference=shared{held}.motion;
if isempty(reference.targets), return; end
k=find(strcmp({reference.targets.name},targetName),1);
if isempty(k)||~reference.clock.phaseReliable, r.status='phase_alignment_unavailable'; return; end
target=reference.targets(k); r.unit=target.unit;
if numel(target.muscles)>1&&~tbench.option(m,'sourceSynchronized',false)
    r.status='EMG_source_clock_unverified'; return;
end
trainX=[]; trainB=[]; trainY=[]; ids={}; sourceCells={}; allNative=true;
levels={'F0','F1','F2','F3','R'}; compressed=struct();
for z=1:numel(levels), compressed.(levels{z})=[]; end
for j=1:numel(shared)
    if j==held||isempty(shared{j}), continue; end
    ref=shared{j}.motion; if isempty(ref.targets), continue; end
    k=find(strcmp({ref.targets.name},targetName),1);
    if isempty(k)||~ref.clock.phaseReliable, continue; end
    t=ref.targets(k); if ~strcmp(t.unit,target.unit)||~strcmp(t.kind,target.kind)||~strcmp(t.frame,target.frame), continue; end
    [X,B,Y]=matrices(packets{j},shared{j},t,m);
    if isempty(B), continue; end
    step=max(1,min(m.waveformStride,floor(t.fs/(4*m.searchBand(2))))); ix=1:step:size(Y,1);
    if isempty(X), allNative=false; else, sourceCells{end+1}=X(ix,:); end %#ok<AGROW>
    if ~isempty(packets{j})&&~isempty(packets{j}.records)
        descriptors=tbench.motion.trialFeatures(packets{j}.records,t.muscles);
        for z=1:numel(levels), key=levels{z}; compressed.(key)=[compressed.(key);repmat(descriptors.(key),numel(ix),1)]; end
    else
        for z=1:numel(levels), key=levels{z}; compressed.(key)=[]; end
    end
    trainB=[trainB;B(ix,:)]; trainY=[trainY;Y(ix,:)]; ids{end+1}=shared{j}.trial; %#ok<AGROW>
end
if allNative&&~isempty(sourceCells), trainX=vertcat(sourceCells{:}); end
if numel(ids)<2, r.status='two_calibration_repetitions_required'; return; end
[X,B,Y]=matrices(packets{held},shared{held},target,m);
good=all(isfinite(Y),2); if sum(good)<m.minimumFitRows, r.status='insufficient_target_support'; return; end
r.testRows=sum(good); r.trainingTrials=strjoin(ids,' '); r.fitRows=sum(all(isfinite(trainY),2));
if isempty(trainB)||isempty(B), r.status='reference_EMG_unavailable'; return; end
predB=zeros(size(Y)); predX=zeros(size(Y)); predFrozen=zeros(size(Y)); native=~isempty(X)&&~isempty(trainX);
trainingResidual=nan(size(trainY));
for d=1:size(Y,2)
    baseline=tbench.motion.fitLinear(trainB,trainY(:,d),m.baselineRidge);
    trainingResidual(:,d)=trainY(:,d)-tbench.motion.predictLinear(baseline,trainB);
    pred=tbench.motion.predictLinear(baseline,B); pred(~isfinite(pred))=0; predB(:,d)=pred;
    if native
        model=tbench.motion.fitLinear(trainX,trainY(:,d),lambda);
        pred=tbench.motion.predictLinear(model,X); pred(~isfinite(pred))=0; predX(:,d)=pred;
        complete=all(isfinite(X),2);
        if any(complete)
            reconstruction=repmat(mean(X(complete,:),1),size(X,1),1);
            pred=tbench.motion.predictLinear(model,reconstruction); pred(~isfinite(pred))=0; predFrozen(:,d)=pred;
        end
    end
end
truth=Y(good,:); errorB=sum((truth-predB(good,:)).^2,2);
r.nullMSE=mean(sum(truth.^2,2)); r.referenceMSE=mean(errorB);
if ~native, r.status='native_waveform_unavailable'; return; end
r.inputCoverage=mean(all(isfinite(X(good,:)),2));
errorX=sum((truth-predX(good,:)).^2,2); r.methodMSE=mean(errorX);
r.gainOverNull=ratio(r.nullMSE-r.methodMSE,r.nullMSE);
r.gainOverReference=ratio(r.referenceMSE-r.methodMSE,r.referenceMSE);
% Common training-fitted reference observation scale; no test-residual fitting.
residual=trainingResidual(isfinite(trainingResidual)); energy=trainY(isfinite(trainY)).^2;
v=0; if ~isempty(residual), v=mean(residual.^2); end
floorVariance=0; if ~isempty(energy), floorVariance=m.modelVarianceFloor*mean(energy); end
v=max([realmin floorVariance v]); r.observationVariance=v;
r.logScoreGainOverReference=mean((errorB-errorX)/(2*v)); r.status='ok';
r.frozenTemporalMeanMSE=mean(sum((truth-predFrozen(good,:)).^2,2));
gain=r.nullMSE-r.methodMSE;
if gain>m.minimumBaselineLoss*r.nullMSE&&gain>0, r.frozenTemporalMeanRetained=(r.nullMSE-r.frozenTemporalMeanMSE)/gain; end
if ~isempty(packets{held}.records)
    descriptors=tbench.motion.trialFeatures(packets{held}.records,target.muscles);
    for z=1:numel(levels)
        key=levels{z}; C=compressed.(key);
        if isempty(C)||size(C,1)~=size(trainY,1), continue; end
        T=repmat(descriptors.(key),size(Y,1),1); prediction=zeros(size(Y));
        for d=1:size(Y,2)
            fit=tbench.motion.fitLinear(C,trainY(:,d),lambda);
            pred=tbench.motion.predictLinear(fit,T); pred(~isfinite(pred))=0; prediction(:,d)=pred;
        end
        loss=mean(sum((truth-prediction(good,:)).^2,2)); label=[key 'waveform'];
        if strcmp(key,'R'), label='RsubjectWaveform'; end
        r.([label 'MSE'])=loss;
        if gain>m.minimumBaselineLoss*r.nullMSE&&gain>0, r.([label 'Retained'])=(r.nullMSE-loss)/gain; end
    end
end
end
function [X,B,Y]=matrices(packet,shared,target,m)
Y=target.waveform; Y(~target.valid,:)=NaN; B=nan(size(Y,1),0); X=B;
if isempty(shared)||~isfield(shared,'baseline')||isempty(shared.baseline), return; end
M=target.muscles; [base,~]=tbench.motion.align(shared.time,shared.baseline(:,M),target,m);
B=quadratures(base);
if isempty(packet)||~all(packet.hasWaveform(M)), return; end
wave=packet.waveform(:,M); wave(~packet.valid(:,M))=NaN;
[wave,~]=tbench.motion.align(packet.time,wave,target,m); X=quadratures(wave);
end
function X=quadratures(y)
X=y;
for k=1:size(y,2), z=tbench.analyticSignal(y(:,k)); X(:,end+1)=imag(z); end %#ok<AGROW>
end
function v=ratio(a,b)
v=NaN; if isfinite(b)&&b>0, v=a/b; end
end
