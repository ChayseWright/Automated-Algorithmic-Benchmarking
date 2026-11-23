function ctx = reference(trial,mvc,cfg)
%REFERENCE Define the power endpoint BEFORE any carrier-rate reduction.
hB=tbench.designFIR(trial.fs,cfg.rawBand,cfg.carrierFIRSeconds);
carrier=tbench.applyFIR(trial.x,hB);
[envelope,t,audit]=tbench.preprocess.rateChange(abs(carrier),trial.fs,cfg.analysisFs, ...
    cfg.reduction.envelopePassHz,cfg.reduction.envelopeStopHz,cfg.reduction);
[energy,~,~]=tbench.preprocess.rateChange(carrier.^2,trial.fs,cfg.analysisFs, ...
    cfg.reduction.envelopePassHz,cfg.reduction.envelopeStopHz,cfg.reduction);
energy(energy<=0)=NaN; valid=isfinite(envelope)&isfinite(energy);
for m=1:size(trial.x,2)
    v=trial.x(isfinite(trial.x(:,m)),m);
    if isempty(v)||max(v)==min(v), valid(:,m)=false; end
end
envelope(~valid)=NaN; energy(~valid)=NaN;
hSlow=tbench.designFIR(cfg.analysisFs,cfg.baselineCutoff,cfg.baselineFIRSeconds);
[slow,ok]=tbench.applyFIR(envelope,hSlow);
ctx=struct('t',t,'fs',cfg.analysisFs,'envelope',envelope,'energy',energy, ...
    'y',envelope-slow,'valid',valid,'dynamicValid',valid&ok,'mvc',mvc, ...
    'rawFs',trial.fs,'raw',trial.x,'rawTime',trial.t,'carrier',carrier, ...
    'carrierFilter',hB,'envelopeFilter',[],'originalDuration',trial.duration);
ctx.observerGuardSeconds=(numel(hB)-1)/(2*trial.fs)+audit.halfSupportSeconds;
ctx.timeOrigin=tbench.option(trial,'timeOrigin',0);
ctx.reductionAudit=audit; ctx.reductionAudit.originalSamples=size(trial.x,1);
ctx.reductionAudit.outputSamples=numel(t);
% Raw carrier candidates retain their own branch; originals stay in input MATs.
factor=max(1,floor(trial.fs/cfg.reduction.minimumCarrierFs+1e-9));
if cfg.reduction.enabled&&factor>1
    rawFs=trial.fs/factor;
    [reduced,tr,rawAudit]=tbench.preprocess.rateChange(trial.x,trial.fs,rawFs, ...
        cfg.rawBand(2),min(rawFs/2*.999, (cfg.rawBand(2)+rawFs/2)/2),cfg.reduction);
else
    reduced=trial.x; tr=trial.t; rawFs=trial.fs;
    rawAudit=struct('inputFs',trial.fs,'outputFs',trial.fs,'ratio',1, ...
        'method','unchanged','halfSupportSeconds',0,'passErrorDB',0,'stopAttenuationDB',Inf);
end
ctx.reducedRaw=reduced; ctx.reducedRawTime=tr; ctx.reducedRawFs=rawFs; ctx.rawReductionAudit=rawAudit;
end
