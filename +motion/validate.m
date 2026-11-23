function trial = validate(record,emgClock,m)
%VALIDATE Signed physical coordinates, explicit units/masks, uniform clock.
if ~isstruct(record)||~isscalar(record)||~all(isfield(record,{'time','targets'}))
    error('tbench:MotionSchema','Each motion trial needs time and targets.');
end
scale=tbench.option(record,'timeScale',m.timeScale);
if ~isnumeric(record.time)||~isreal(record.time)||~isvector(record.time)|| ...
        ~isnumeric(scale)||~isscalar(scale)||~isfinite(scale)||scale<=0
    error('tbench:MotionTime','Time must be a real vector with a positive finite timeScale.');
end
t=double(record.time(:))*scale;
if numel(t)<3||any(~isfinite(t))||any(diff(t)<=0), error('tbench:MotionTime','Motion time must be finite and strictly increasing.'); end
sync=tbench.option(record,'sync',struct()); c=tbench.motion.clock(sync,t,emgClock.timeOrigin,m);
original=t;
if c.aligned, t=c.offsetSeconds+c.rateRatio*t-emgClock.timeOrigin; else, t=t-t(1); end
dt=median(diff(t)); fs=1/dt;
if fs<m.minimumFs*(1-1e-6), error('tbench:MotionRate','Motion Fs is below the declared minimum for the target band.'); end
ix=round((t-t(1))/dt)+1;
if any(abs(t-(t(1)+(ix-1)*dt))>m.timeTolerance*dt)||numel(unique(ix))~=numel(ix)
    error('tbench:MotionJitter','Motion timestamps exceed the uniform-clock tolerance.');
end
if ix(end)>m.maxGapExpansion*numel(t), error('tbench:MotionGaps','Motion gap expansion exceeds the declared bound.'); end
targets=record.targets;
if ~isstruct(targets)||isempty(targets)||~all(isfield(targets,{'name','kind','values','unit'}))
    error('tbench:MotionTargets','Provide named signed-coordinate targets with kind, values and unit.');
end
for j=1:numel(targets)
    if ~ischar(targets(j).name)||isempty(targets(j).name)||size(targets(j).name,1)~=1
        error('tbench:MotionTargets','Each target needs a nonempty character-row name.');
    end
end
if numel(unique({targets.name}))~=numel(targets), error('tbench:MotionTargets','Motion target names must be unique.'); end
normalized=struct([]);
for j=1:numel(targets)
    s=targets(j);
    if ~all(isfield(s,{'name','kind','values','unit'})), error('tbench:MotionTarget','Targets need name, kind, values and unit.'); end
    if ~isnumeric(s.values)||~isreal(s.values), error('tbench:MotionShape','Target values must be real numeric coordinates.'); end
    v=double(s.values);
    if ~isreal(v)||size(v,1)~=numel(t)||isempty(v)||size(v,2)>3
        error('tbench:MotionShape','A target needs N-by-1..3 signed coordinate values.');
    end
    kind=char(s.kind); unit=char(s.unit); factor=1;
    if strcmp(kind,'translation')
        if strcmp(unit,'m'), factor=1000; elseif ~strcmp(unit,'mm'), error('tbench:MotionUnit','Translation units must be m or mm.'); end
        unit='mm';
    elseif any(strcmp(kind,{'rotation','joint_angle'}))
        if strcmp(unit,'deg'), factor=pi/180; elseif ~strcmp(unit,'rad'), error('tbench:MotionUnit','Angular units must be rad or deg.'); end
        unit='rad';
    else, error('tbench:MotionKind','Use translation, rotation or joint_angle.');
    end
    quality=flags(tbench.option(s,'valid',true(size(v))),'valid');
    if isvector(quality)&&numel(quality)==numel(t), quality=repmat(logical(quality(:)),1,size(v,2)); end
    if ~isequal(size(quality),size(v)), error('tbench:MotionQuality','Target valid must be N-by-1 or N-by-D.'); end
    good=logical(quality)&isfinite(v); v=v*factor; v(~good)=NaN;
    if strcmp(kind,'joint_angle')
        for d=1:size(v,2)
            runs=tbench.validRuns(isfinite(v(:,d)));
            for z=1:size(runs,1), a=runs(z,1):runs(z,2); v(a,d)=unwrap(v(a,d)); end
        end
    end
    values=nan(ix(end),size(v,2)); values(ix,:)=v;
    conf=flags(tbench.option(s,'confounded',false(numel(t),1)),'confounded');
    if numel(conf)~=numel(t), error('tbench:MotionConfound','confounded must have one flag per native sample.'); end
    ambiguous=true(ix(end),1); ambiguous(ix)=logical(conf(:));
    sd=tbench.option(s,'noiseSD',[]);
    if ~isnumeric(sd)||~isreal(sd), error('tbench:MotionNoise','noiseSD must be numeric and real.'); end
    if ~isempty(sd), sd=double(sd(:)')*factor; end
    if ~isempty(sd)&&~(numel(sd)==1||numel(sd)==size(v,2)), error('tbench:MotionNoise','noiseSD is scalar or one value per coordinate.'); end
    if any(~isfinite(sd)|sd<=0), error('tbench:MotionNoise','Declared noiseSD must be positive and finite.'); end
    limit=tbench.option(s,'limitRMS',NaN);
    if ~isnumeric(limit)||~isreal(limit), error('tbench:MotionLimit','limitRMS must be numeric and real.'); end
    limit=limit*factor;
    if ~isscalar(limit)||(~isnan(limit)&&(~isfinite(limit)||limit<=0)), error('tbench:MotionLimit','limitRMS is a positive physical RMS amplitude or NaN.'); end
    muscles=tbench.option(s,'muscles',1:15);
    if isempty(muscles)||any(muscles<1|muscles>15|muscles~=floor(muscles))||numel(unique(muscles))~=numel(muscles)
        error('tbench:MotionMuscles','Target muscles must be unique acquisition indices 1..15.');
    end
    row=struct('name',char(s.name),'kind',kind,'unit',unit,'frame',char(tbench.option(s,'frame','unspecified')), ...
        'values',values,'confounded',ambiguous,'noiseSD',sd,'limitRMS',limit, ...
        'calibrationVerified',scalarFlag(tbench.option(s,'calibrationVerified',false),'calibrationVerified'), ...
        'calibrationPolicyHash',char(tbench.option(s,'calibrationPolicyHash','')), ...
        'reviewed',scalarFlag(tbench.option(s,'reviewed',tbench.option(record,'reviewed',false)),'reviewed'), ...
        'muscles',muscles(:)','qualitySupplied',isfield(s,'valid'));
    if isempty(normalized), normalized=row; else, normalized(end+1)=row; end %#ok<AGROW>
end
trial=struct('time',t(1)+(0:ix(end)-1)'*dt,'fs',fs,'targets',normalized,'clock',c, ...
    'originalStartTime',original(1),'missingSamples',ix(end)-numel(t));
end
function v=flags(v,name)
if ~(isnumeric(v)||islogical(v))||~isreal(v)||any(~isfinite(v(:)))||any(v(:)~=0&v(:)~=1)
    error('tbench:MotionQuality','%s must contain finite logical or numeric 0/1 flags.',name);
end
v=logical(v);
end
function v=scalarFlag(v,name)
v=flags(v,name); if ~isscalar(v), error('tbench:MotionQuality','%s must be scalar.',name); end
end
