function out = run(prepared,cfg,pipeline)
%RUN Analysis may see one muscle, never clinical labels or other muscles.
local=cfg; p=pipeline.preprocessParams; a=pipeline.analysisParams;
local.spectralWindow=p.windowSeconds; local.localWindow=p.localWindowSeconds;
local.localHop=p.localHop; local.spectralOverlap=p.spectralOverlap;
if strcmp(pipeline.inputKind,'carrier_time'), local.rawBand(1)=p.rawHighpassHz; end
if isfield(a,'threshold'), local.cfar.ratio=a.threshold; end
if isfield(a,'minCycles'), local.local.minCycles=a.minCycles; end
if isfield(a,'switchPenalty'), local.local.switchPenalty=a.switchPenalty; end
if isfield(a,'frequencyPenalty'), local.local.frequencyPenalty=a.frequencyPenalty; end
if isfield(a,'NW'), local.local.NW=a.NW; end
if isfield(a,'nTapers'), local.local.nTapers=a.nTapers; end
if isfield(a,'cfarGuardHz'), local.cfar.guardHz=a.cfarGuardHz; end
if isfield(a,'cfarReferenceHz'), local.cfar.referenceHz=a.cfarReferenceHz; end
if isfield(a,'cfarQuantile'), local.cfar.quantile=a.cfarQuantile; end
out=tbench.methodProposed(prepared,local,a,pipeline.analysisId);
if strcmp(out.status,'ok')
    native=out; reference=prepared;
    reference.envelope=prepared.referenceEnvelope; reference.energy=prepared.referenceEnergy;
    reference.y=prepared.referenceY; reference.valid=prepared.referenceValid;
    reference.dynamicValid=prepared.referenceDynamicValid;
    if out.capabilities.dynamic
        out=tbench.proposalOutput(reference,local,out.frequencyTrack,out.active,out.valid,native.audit);
        out.onProbability=native.onProbability; out.capabilities.probabilistic=native.capabilities.probabilistic;
        out.audit.powerEstimator='fixed_reference_harmonic_refit';
    else
        f=tbench.option(native.audit,'proposalFrequency',native.frequency);
        % The same fixed Welch power meter avoids global line-fit attenuation
        % from an otherwise small frequency-bin error. Detector statistics differ.
        out=tbench.analysis.welchPower(reference,cfg,f,native.detected,native.audit);
        out.spectrum=native.spectrum; out.spectrumFrequency=native.spectrumFrequency;
    end
    out.audit.nativeDetectorPower=native.powerMean;
    out.audit.nativeDetectorFrequency=native.frequency;
    out.audit.nativePowerEstimator='method_internal_before_common_reporting_adapter';
end
out=tbench.validateOutput(out,prepared);
out.audit.preprocessing=prepared.preprocessAudit;
out.audit.pipeline=pipeline.id; out.audit.candidate=pipeline.candidateId;
end
