function rows = calibrate(candidates,cfg,seedOffset)
%CALIBRATE Independent raw-Gaussian AM challenges in the fixed reporting domain.
if nargin<3, seedOffset=0; end
if isempty(candidates), rows=struct([]); return; end
state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
rng(cfg.seed+100+seedOffset,'twister'); rows=struct([]);
fs=cfg.calibration.fs; T=cfg.calibration.duration; t=(0:round(T*fs)-1)'/fs;
hB=tbench.designFIR(fs,cfg.rawBand,cfg.carrierFIRSeconds);
freq=cfg.calibration.frequencies; depths=cfg.calibration.depths;
ppKeys=cell(1,numel(candidates)); detectorKeys=ppKeys;
for a=1:numel(candidates)
    ppKeys{a}=tbench.study.preparationKey(candidates(a)); p=candidates(a).preprocessParams;
    detectorKeys{a}=sprintf('%.17g_%.17g_%.17g',p.rawHighpassHz,p.envelopePassHz,p.envelopeStopHz);
end
[~,~,ppIndex]=unique(ppKeys); [~,~,detectorIndex]=unique(detectorKeys);
lastPP=zeros(1,max(ppIndex)); lastDetector=zeros(1,max(detectorIndex));
for a=1:numel(candidates), lastPP(ppIndex(a))=a; lastDetector(detectorIndex(a))=a; end
for k=1:cfg.calibration.nCases
    f=freq(mod(k-1,numel(freq))+1);
    depth=depths(mod(floor((k-1)/numel(freq)),numel(depths))+1);
    patterns=tbench.option(cfg.calibration,'patterns',{'continuous','middle'});
    pattern=patterns{mod(floor((k-1)/(numel(freq)*numel(depths))),numel(patterns))+1};
    intermittent=~strcmp(pattern,'continuous'); gate=ones(size(t)); duration=max(3/f,.5);
    starts=[];
    if intermittent
        switch pattern
            case 'early', starts=.20*T;
            case 'middle', starts=.45*T;
            case 'late', starts=.80*T;
            case 'multiple', starts=[.20 .50 .80]*T;
            otherwise, error('tbench:CalibrationPattern','Unknown episode challenge.');
        end
        gate(:)=0;
        for start=starts, gate=gate+double(t>=start&t<min(T,start+duration)); end
        gate=min(gate,1);
    end
    drift=0;
    if mod(k,4)==0&&f>cfg.searchBand(1)+.25&&f<cfg.searchBand(2)-.25, drift=.5; end
    phase=2*pi*(f*t+drift*(t.^2/(2*T)-t/2)); modulation=1+depth*gate.*cos(phase);
    noise=.3; raw=modulation.*randn(size(t))+noise*randn(size(t));
    transient=mod(k,7)==0;
    if transient, raw=raw+10*exp(-.5*((t-.7*T)/.015).^2); end
    local=tbench.motion.inferenceConfig(cfg); local.timeScale=1; local.nominalFs=[]; local.knownRails=[];
    gap=mod(k,11)==0; clipped=mod(k,13)==0;
    if gap, raw(t>=.7*T&t<.7*T+.2)=NaN; end
    if clipped
        raw=max(-.8,min(.8,raw)); local.knownRails=[-.8 .8];
    end
    if clipped, raw(raw<=-.8|raw>=.8)=NaN; end
    trial=struct('x',raw,'t',t,'fs',fs,'duration',numel(t)/fs,'timeOrigin',0);
    ctx=tbench.preprocess.reference(trial,1,local);
    variance=conv(modulation.^2+noise^2,hB.^2,'same');
    expectedRaw=sqrt(2/pi)*sqrt(max(variance,0));
    [truthEnvelope,~,~]=tbench.preprocess.rateChange(expectedRaw,fs,ctx.fs, ...
        cfg.reduction.envelopePassHz,cfg.reduction.envelopeStopHz,cfg.reduction);
    truthSlow=tbench.applyFIR(truthEnvelope,tbench.designFIR(ctx.fs,cfg.baselineCutoff,cfg.baselineFIRSeconds));
    truthY=truthEnvelope-truthSlow;
    trueGate=interp1(t,gate,ctx.t,'nearest',0)>0;
    phi=2*pi*(f*ctx.t+drift*(ctx.t.^2/(2*T)-ctx.t/2));
    B=[ones(size(phi)) double(trueGate) double(trueGate).*cos(phi) double(trueGate).*sin(phi)];
    valid=ctx.dynamicValid(:,1)&isfinite(truthY); beta=pinv(B(valid,:))*truthY(valid);
    truthFund=B(:,3:4)*beta(3:4);
    if depth==0, truthFund(:)=0; end
    c=tbench.study.channel(ctx,1);
    prepared=cell(size(lastPP)); detectors=cell(size(lastDetector));
    for a=1:numel(candidates)
        begin=tic;
        try
            q=ppIndex(a); d=detectorIndex(a);
            if isempty(prepared{q})
                p=candidates(a).preprocessParams;
                if isempty(detectors{d}), detectors{d}=tbench.preprocess.detector(c,local,p); end
                input=c; input.cachedDetector=detectors{d};
                input.cachedDetectorParameters=[p.rawHighpassHz p.envelopePassHz p.envelopeStopHz];
                prepared{q}=candidates(a).preprocess(input,local,candidates(a));
            end
            out=candidates(a).analyze(prepared{q},local,candidates(a));
        catch err
            if cfg.failFast, rethrow(err); end
            out=tbench.emptyOutput(c); out.status='method_error'; out.message=err.message;
        end
        ok=out.valid&valid; truthPower=NaN; if any(ok), truthPower=mean(truthFund(ok).^2); end
        estimated=out.powerMean;
        error=NaN; if depth>0&&isfinite(out.frequency), error=abs(out.frequency-f); end
        relative=NaN; if truthPower>0&&isfinite(estimated), relative=abs(estimated-truthPower)/truthPower; end
        tp=NaN; fp=NaN; fn=NaN; iou=NaN; onset=NaN; offset=NaN; trackRMSE=NaN; activeTruthPower=NaN; activeEstimatePower=NaN;
        if out.capabilities.dynamic&&~isempty(out.active)
            actual=trueGate&ok&depth>0; found=logical(out.active)&ok;
            tp=sum(actual&found); fp=sum(~actual&found&ok); fn=sum(actual&~found);
            if tp+fp+fn>0, iou=tp/(tp+fp+fn); end
            if any(actual), activeTruthPower=mean(truthFund(actual).^2); activeEstimatePower=mean(out.fundamental(actual).^2); end
            common=actual&found&isfinite(out.frequencyTrack);
            if any(common), trackRMSE=sqrt(mean((out.frequencyTrack(common)-(f+drift*(ctx.t(common)/T-.5))).^2)); end
            if any(actual)&&any(found)
                onset=ctx.t(find(found,1))-ctx.t(find(actual,1));
                offset=ctx.t(find(found,1,'last'))-ctx.t(find(actual,1,'last'));
            end
        end
        row=struct('method',candidates(a).analysisId,'pipelineId',candidates(a).id, ...
            'candidateId',candidates(a).candidateId,'candidateHash',candidates(a).hash, ...
            'caseIndex',k,'seedOffset',seedOffset,'status',out.status, ...
            'truthScope','artifact_free_expected_reporting_envelope_projection', ...
            'frequency',f,'depth',depth,'intermittent',intermittent,'pattern',pattern,'transient',transient,'gap',gap,'clipped',clipped, ...
            'detected',out.detected,'truthPower',truthPower,'estimatedPower',estimated, ...
            'relativePowerError',relative,'frequencyAbsError',error,'TPsamples',tp,'FPsamples',fp,'FNsamples',fn, ...
            'episodeIoU',iou,'onsetSignedError',onset,'offsetSignedError',offset, ...
            'trackRMSEConditional',trackRMSE,'truthActivePower',activeTruthPower, ...
            'estimatedOnTruthActivePower',activeEstimatePower, ...
            'validFraction',mean(ok),'runtimeSeconds',toc(begin));
        if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
        if a==lastPP(ppIndex(a)), prepared{ppIndex(a)}=[]; end
        if a==lastDetector(detectorIndex(a)), detectors{detectorIndex(a)}=[]; end
    end
end
end
