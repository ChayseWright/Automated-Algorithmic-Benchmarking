function value = policyHash(m)
%POLICYHASH Calibration applies to this reference policy, not readout tuning.
names={'maxFs','minimumFs','slowCutoffHz','slowSeconds','searchBand','witnessBand', ...
    'windowSeconds','shortWindowSeconds','hopSeconds','NW','nTapers','frequencyStepHz', ...
    'timeScale','timeTolerance','maxGapExpansion','requireKineticReview', ...
    'minimumCycles','minimumEpisodeSeconds','evidenceRatio','frequencyAgreementHz','confidenceZ','reduction','oscillator'};
s=struct(); for k=1:numel(names), s.(names{k})=m.(names{k}); end
s.schemaVersion=2;
functions={'tbench.motion.reference','tbench.motion.validate','tbench.motion.oscillator', ...
    'tbench.motion.harmonicFit','tbench.motion.sample','tbench.motion.clock','tbench.tapers', ...
    'tbench.designFIR','tbench.applyFIR','tbench.preprocess.rateChange','tbench.preprocess.designAA'};
for k=1:numel(functions)
    path=which(functions{k});
    if isempty(path), error('tbench:MotionReferenceCode','Reference implementation is incomplete: %s.',functions{k}); end
    s.implementation.(sprintf('source_%02d',k))=tbench.study.fileHash(path);
end
value=tbench.study.hash(s);
end
