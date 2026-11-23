function validateConfig(m,cfg)
%VALIDATECONFIG Reject impossible motion policies before reading recordings.
flags={'enabled','required','requireKineticReview','selectionTieBreak'};
for k=1:numel(flags), boolean(m.(flags{k}),flags{k}); end
if ~m.enabled, return; end
positive={'maxFs','minimumFs','timeScale','timeTolerance','maxGapExpansion','slowCutoffHz', ...
    'slowSeconds','windowSeconds','shortWindowSeconds','hopSeconds','NW','frequencyStepHz', ...
    'minimumCycles','minimumEpisodeSeconds','evidenceRatio','frequencyAgreementHz', ...
    'confidenceZ','minimumBaselineLoss','modelVarianceFloor','baselineRidge'};
for k=1:numel(positive)
    v=m.(positive{k});
    if ~isnumeric(v)||~isreal(v)||~isscalar(v)||~isfinite(v)||v<=0
        error('tbench:MotionConfig','motion.%s must be a positive finite scalar.',positive{k});
    end
end
fractions={'minimumReferenceCoverage','minimumSourceCoverage'};
for k=1:numel(fractions)
    v=m.(fractions{k});
    if ~isnumeric(v)||~isscalar(v)||~isfinite(v)||v<=0||v>1
        error('tbench:MotionConfig','motion.%s must be in (0,1].',fractions{k});
    end
end
integers={'nTapers','minimumActiveMuscles','waveformStride','minimumFitRows','minimumMotionSubjects'};
for k=1:numel(integers)
    v=m.(integers{k});
    if ~isnumeric(v)||~isscalar(v)||~isfinite(v)||v<1||v~=floor(v)
        error('tbench:MotionConfig','motion.%s must be a positive integer.',integers{k});
    end
end
for key={'searchBand','witnessBand'}
    v=m.(key{1});
    if ~isnumeric(v)||numel(v)~=2||any(~isfinite(v))||v(1)<=0||v(2)<=v(1)
        error('tbench:MotionConfig','Motion bands need increasing positive endpoints.');
    end
end
if ~isequal(m.searchBand(:)',cfg.searchBand(:)')
    error('tbench:MotionConfig','Primary EMG and motion search bands must match for this comparison.');
end
if m.maxFs<m.minimumFs||.45*m.minimumFs<=1.1*m.searchBand(2)|| ...
        m.slowCutoffHz>=m.searchBand(1)||m.witnessBand(1)>m.searchBand(1)|| ...
        m.witnessBand(2)<m.searchBand(2)||m.shortWindowSeconds>m.windowSeconds|| ...
        m.hopSeconds>m.shortWindowSeconds||m.nTapers>floor(2*m.NW-1)||m.minimumActiveMuscles>15
    error('tbench:MotionConfig','Incompatible motion rate, band, window or taper policy.');
end
if ~iscellstr(m.primaryTargets)||numel(unique(m.primaryTargets))~=numel(m.primaryTargets)|| ...
        ~any(strcmp(m.clinicalTargetKind,{'translation','rotation','joint_angle'}))
    error('tbench:MotionConfig','Use unique target names and one physical coordinate kind for the clinical comparator.');
end
for key={'maxTimingErrorSeconds','maxPhaseErrorDegrees'}
    v=m.(key{1});
    if ~isnumeric(v)||~isscalar(v)||~isfinite(v)||v<0
        error('tbench:MotionConfig','Timing/phase tolerances must be finite and nonnegative.');
    end
end
if isempty(m.readoutRidges)||~isnumeric(m.readoutRidges)||any(~isfinite(m.readoutRidges)|m.readoutRidges<=0)|| ...
        numel(unique(m.readoutRidges))~=numel(m.readoutRidges)
    error('tbench:MotionConfig','Use unique positive ridge values for mechanical readouts.');
end
c=m.coherence; boolean(c.enabled,'coherence.enabled'); boolean(m.calibration.enabled,'calibration.enabled');
if c.windowSeconds<=0||c.NW<1||c.nTapers<1||c.nTapers~=floor(c.nTapers)||c.nTapers>floor(2*c.NW-1)|| ...
        c.minimumSegments<2||c.minimumSegments~=floor(c.minimumSegments)||c.replicatesPerInput<1|| ...
        c.maximumCondition<=1||c.shrinkage<=0||c.shrinkage>1|| ...
        c.minimumMultipleCoherence<0||c.minimumMultipleCoherence>1||c.surrogates<0||c.surrogates~=floor(c.surrogates)
    error('tbench:MotionConfig','Invalid coherence capacity/identifiability settings.');
end
numeric=[c.windowSeconds c.NW c.nTapers c.minimumSegments c.replicatesPerInput ...
    c.maximumCondition c.shrinkage c.minimumMultipleCoherence c.surrogates c.seed];
if any(~isfinite(numeric)), error('tbench:MotionConfig','Coherence settings must be finite.'); end
o=m.oscillator;
if any(~isfinite([o.damping o.processFraction o.stayOn o.stayOff o.onThreshold]))|| ...
        o.damping<=0||o.damping>1||o.processFraction<=0|| ...
        o.stayOn<=0||o.stayOn>=1||o.stayOff<=0||o.stayOff>=1||o.onThreshold<=.5||o.onThreshold>=1
    error('tbench:MotionConfig','Invalid switching-oscillator reference settings.');
end
if m.reduction.stopAttenuationDB<60||m.reduction.passRippleDB<=0|| ...
        m.reduction.sincHalfLength<8||m.reduction.sincHalfLength~=floor(m.reduction.sincHalfLength)|| ...
        any(~isfinite([m.reduction.stopAttenuationDB m.reduction.passRippleDB m.reduction.sincHalfLength]))
    error('tbench:MotionConfig','Invalid motion anti-alias/resampling specification.');
end
if m.calibration.enabled
    a=m.calibration;
    if a.fs<m.minimumFs||a.duration<2*m.windowSeconds||a.noiseSD<=0|| ...
            any(~isfinite([a.fs a.duration a.noiseSD a.frequencies a.amplitudesRMS]))|| ...
            isempty(a.frequencies)||isempty(a.amplitudesRMS)||any(a.amplitudesRMS<0)|| ...
            any(a.frequencies<m.searchBand(1)|a.frequencies>m.searchBand(2))
        error('tbench:MotionConfig','Invalid synthetic motion challenge settings.');
    end
end
end
function boolean(v,name)
if ~(islogical(v)||isnumeric(v))||~isscalar(v)||~isfinite(v)||~any(v==[0 1])
    error('tbench:MotionConfig','motion.%s must be a scalar boolean.',name);
end
end
