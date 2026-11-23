function c = clock(sync,time,emgOrigin,m)
%CLOCK Hardware-event/declared affine alignment; never EMG correlation.
c=struct('aligned',false,'phaseReliable',false,'timingReliable',false, ...
    'offsetSeconds',NaN,'rateRatio',NaN,'uncertaintySeconds',Inf, ...
    'maximumPhaseErrorDegrees',Inf,'method','unverified_independent_trial_clock', ...
    'timeOrigin',emgOrigin,'eventCount',0,'eventResidualSeconds',NaN);
if nargin<1||isempty(sync), return; end
if ~isstruct(sync)||~isscalar(sync), error('tbench:MotionSync','sync must be a scalar structure.'); end
unc=tbench.option(sync,'uncertaintySeconds',Inf);
if ~isnumeric(unc)||~isreal(unc)||~isscalar(unc)||unc<0||isnan(unc), error('tbench:MotionSync','Timing uncertainty must be nonnegative.'); end
verified=tbench.option(sync,'verified',false);
if ~(isnumeric(verified)||islogical(verified))||~isscalar(verified)||~isfinite(verified)||~any(verified==[0 1])
    error('tbench:MotionSync','sync.verified must be a finite scalar boolean.');
end
a=tbench.option(sync,'motionEvents',[]); b=tbench.option(sync,'emgEvents',[]);
if ~isnumeric(a)||~isnumeric(b)||~isreal(a)||~isreal(b)
    error('tbench:MotionSync','Hardware event times must be real numeric seconds.');
end
if ~isempty(a)||~isempty(b)
    a=a(:); b=b(:);
    if numel(a)~=numel(b)||any(~isfinite([a;b]))||isempty(a)
        error('tbench:MotionSync','Provide equally sized finite event arrays in seconds.');
    end
    c.eventCount=numel(a);
    if numel(a)>=2&&max(a)>min(a)
        center=mean(a); X=[ones(size(a)) a-center]; fit=X\b;
        c.offsetSeconds=fit(1)-fit(2)*center; c.rateRatio=fit(2);
        c.eventResidualSeconds=max(abs(X*fit-b)); c.method='hardware_events_affine';
        unc=unc+c.eventResidualSeconds;
    else
        ppm=tbench.option(sync,'clockPPMbound',Inf);
        if ~isnumeric(ppm)||~isscalar(ppm)||isnan(ppm)||ppm<0
            error('tbench:MotionSync','clockPPMbound must be nonnegative.');
        end
        c.offsetSeconds=b(1)-a(1); c.rateRatio=1;
        unc=unc+max(abs(time-a(1)))*ppm*1e-6; c.method='single_event_with_drift_bound';
    end
elseif verified
    c.offsetSeconds=tbench.option(sync,'offsetSeconds',NaN);
    c.rateRatio=tbench.option(sync,'rateRatio',NaN); c.method='declared_verified_affine';
end
if ~isnumeric(c.offsetSeconds)||~isnumeric(c.rateRatio)||~isscalar(c.offsetSeconds)||~isscalar(c.rateRatio)||~isreal([c.offsetSeconds c.rateRatio])
    error('tbench:MotionSync','Affine clock offset and rate ratio must be real scalars.');
end
if isfinite(c.offsetSeconds)&&isfinite(c.rateRatio)&&c.rateRatio>0&&isfinite(unc)
    c.aligned=true; c.uncertaintySeconds=unc;
    c.maximumPhaseErrorDegrees=360*m.searchBand(2)*unc;
    c.timingReliable=unc<=m.maxTimingErrorSeconds;
    c.phaseReliable=c.maximumPhaseErrorDegrees<=m.maxPhaseErrorDegrees;
end
end
