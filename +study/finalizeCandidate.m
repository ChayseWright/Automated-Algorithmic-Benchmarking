function p = finalizeCandidate(p,cfg)
%FINALIZECANDIDATE Validate theory constraints and assign a content identity.
a=p.analysisParams; pp=p.preprocessParams;
if pp.rawHighpassHz<0||pp.rawHighpassHz>=cfg.rawBand(2)||pp.envelopePassHz<cfg.harmonics*cfg.searchBand(2)|| ...
        pp.envelopeStopHz>=cfg.analysisFs/2||pp.envelopePassHz>=pp.envelopeStopHz
    error('tbench:CandidateBandwidth','Invalid detector carrier/envelope passband.');
end
if isfield(a,'nTapers')&&isfield(a,'NW')
    % Repair an inadmissible taper pair to its smallest valid bandwidth.
    a.NW=max(a.NW,(a.nTapers+1)/2);
end
if strcmp(p.inputKind,'complex_tf')&&isfield(a,'cycles')
    gamma=tbench.option(a,'morseGamma',3); a.morseGamma=gamma; a.morseBeta=a.cycles^2/gamma;
end
if isfield(a,'minCycles')&&any(strcmp(p.inputKind,{'local_spectrum','envelope_time'}))&& ...
        any(strcmp(p.analysisFamily,{'legacy','representation','auxiliary'}))
    pp.localWindowSeconds=max(pp.localWindowSeconds,a.minCycles/cfg.searchBand(1));
end
if pp.windowSeconds<=0||pp.localWindowSeconds<=0||pp.localHop<=0||pp.spectralOverlap<0||pp.spectralOverlap>=1
    error('tbench:CandidateSupport','Window/hop/overlap parameters are inadmissible.');
end
if isfield(a,'threshold')&&(~isfinite(a.threshold)||a.threshold<=0), error('tbench:Threshold','Evidence threshold must be positive.'); end
p.analysisParams=a; p.preprocessParams=pp;
p.hash=tbench.study.hash(struct('pipelineIdentity',p.id,'implementation',p.implementationHash,'preprocess',p.preprocessId, ...
    'preprocessOwner',p.preprocessAnalysisId,'pp',pp, ...
    'analysis',p.analysisId,'ap',a));
p.candidateId=['candidate_' p.hash(1:12)];
end
