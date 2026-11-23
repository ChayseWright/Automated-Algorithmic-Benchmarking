function ref = reference(record,emgClock,m,task)
%REFERENCE Independent kinematic target; input contains no EMG or diagnosis.
trial=tbench.motion.validate(record,emgClock,m); targets=struct([]); policy=tbench.motion.policyHash(m);
fsOut=min(m.maxFs,trial.fs); stop=.45*fsOut;
pass=max(1.10*m.searchBand(2),min(30,.36*fsOut));
if pass>=stop, error('tbench:MotionBandwidth','Motion rate cannot preserve the primary band with the declared anti-alias margin.'); end
for j=1:numel(trial.targets)
    s=trial.targets(j); origin=trial.time(1);
    [q,t,audit]=tbench.preprocess.rateChange(s.values,trial.fs,fsOut,pass,stop,m.reduction); t=t+origin;
    slow=tbench.applyFIR(q,tbench.designFIR(fsOut,m.slowCutoffHz,m.slowSeconds)); y=q-slow;
    nominalSeconds=tbench.option(emgClock,'duration',numel(t)/fsOut);
    if ~isnumeric(nominalSeconds)||~isscalar(nominalSeconds)||~isfinite(nominalSeconds)||nominalSeconds<=0
        error('tbench:MotionClock','Nominal trial duration must be positive finite seconds.');
    end
    if trial.clock.aligned&&isfinite(nominalSeconds)&&nominalSeconds>0
        use=t>=0&t<nominalSeconds;
        q=q(use,:); slow=slow(use,:); y=y(use,:); t=t(use);
    end
    if numel(t)<round(m.shortWindowSeconds*fsOut)
        error('tbench:MotionOverlap','Insufficient motion support in the nominal EMG trial interval.');
    end
    valid=all(isfinite(y),2);
    y(~valid,:)=NaN; % All vector axes use exactly the same physical support.
    conf=interp1(trial.time,double(s.confounded),t,'nearest',1)>.5;
    witnessBand=[max(.5,m.witnessBand(1)) min(m.witnessBand(2),pass-.25)];
    witness=tbench.applyFIR(y,tbench.designFIR(fsOut,witnessBand,2));
    [P,f]=vectorSpectrum(y,fsOut,min(4,numel(t)/fsOut),.5,m.NW,m.nTapers,m.frequencyStepHz);
    integrationBand=[max(.5,m.searchBand(1)-.5) min(pass,m.searchBand(2)+.5)];
    band=f>=integrationBand(1)&f<=integrationBand(2); meanPower=NaN;
    if sum(band)>=2&&any(isfinite(P(band))), meanPower=trapz(f(band),P(band)); end
    secondary=tbench.motion.oscillator(y,valid,fsOut,m);
    N=numel(t); wPower=nan(N,1); wLow=wPower; wHigh=wPower; track=wPower; ratio=wPower;
    supported=false(N,1); frames=struct([]);
    windows=[m.windowSeconds m.shortWindowSeconds];
    centers=(1:max(1,round(m.hopSeconds*fsOut)):N)';
    for ii=1:numel(centers)
        center=centers(ii); choices=struct([]);
        for seconds=windows
            n=max(20,round(seconds*fsOut)); a=center-floor(n/2); b=a+n-1;
            if a<1||b>N||~all(valid(a:b)), continue; end
            [S,g]=vectorSpectrum(y(a:b,:),fsOut,n/fsOut,0,m.NW,m.nTapers,m.frequencyStepHz);
            eligible=find(g>=m.searchBand(1)&g<=m.searchBand(2)&isfinite(S));
            if isempty(eligible), continue; end
            bg=g>=m.witnessBand(1)&g<=min(m.witnessBand(2),pass); noise=nan(size(S));
            for k=eligible(:)'
                use=bg&abs(g-g(k))>.75&isfinite(S);
                z=S(use); if ~isempty(z), noise(k)=median(z); end
            end
            evidence=S./max(noise,realmin); [value,index]=max(evidence(eligible)); index=eligible(index); frequency=g(index);
            if ~isfinite(value)||seconds*frequency<m.minimumCycles, continue; end
            if index>1&&index<numel(g)&&all(S(index-1:index+1)>0)
                v=log(S(index-1:index+1)); denominator=v(1)-2*v(2)+v(3);
                if abs(denominator)>1e-12
                    delta=max(-.5,min(.5,.5*(v(1)-v(3))/denominator)); frequency=frequency+delta*(g(2)-g(1));
                end
            end
            frequency=max(m.searchBand(1),min(m.searchBand(2),frequency));
            fit=tbench.motion.harmonicFit(y(a:b,:),t(a:b),frequency,fsOut);
            row=struct('time',t(center),'frequency',frequency,'powerSigned',fit.powerSigned, ...
                'powerRaw',fit.powerRaw,'powerSD',fit.powerSD,'evidence',value,'seconds',seconds);
            if isempty(choices), choices=row; else, choices(end+1)=row; end %#ok<AGROW>
        end
        if isempty(choices), continue; end
        % The long witness is preferred; a short witness is used when stronger.
        [~,pick]=max([choices.evidence]); row=choices(pick);
        if isempty(frames), frames=row; else, frames(end+1)=row; end %#ok<AGROW>
    end
    if numel(frames)>=2
        a=[frames.time]'; nearest=interp1(a,(1:numel(a))',t,'nearest',NaN); good=isfinite(nearest); ix=nearest(good);
        p=[frames.powerSigned]; sd=[frames.powerSD]; ff=[frames.frequency]; e=[frames.evidence];
        wPower(good)=p(ix); wLow(good)=p(ix)-m.confidenceZ*sd(ix); wHigh(good)=p(ix)+m.confidenceZ*sd(ix);
        track(good)=ff(ix); ratio(good)=e(ix);
        supported(good)=abs(t(good)-a(ix))<=m.hopSeconds/2+1/fsOut;
    end
    supported=supported&valid;
    wPower(~supported)=NaN; wLow(~supported)=NaN; wHigh(~supported)=NaN;
    track(~supported)=NaN; ratio(~supported)=NaN;
    reviewOK=task(1)~='K'||~m.requireKineticReview||s.reviewed;
    coherent=abs(track-secondary.frequency)<=m.frequencyAgreementHz;
    positive=supported&~conf&reviewOK&ratio>=m.evidenceRatio&wLow>0& ...
        secondary.onScore>=m.oscillator.onThreshold&coherent;
    state=zeros(N,1,'int8'); state(positive)=1;
    calibrated=s.calibrationVerified&&strcmp(s.calibrationPolicyHash,policy)&& ...
        ~isempty(s.noiseSD)&&isfinite(s.limitRMS);
    negative=supported&~conf&reviewOK&calibrated&wHigh<s.limitRMS^2&secondary.onScore<.5;
    state(negative)=-1;
    positive=minimumEvents(state==1,t,track,fsOut,m); state(state==1&~positive)=0;
    multitaperState=zeros(N,1,'int8'); oscillatorState=multitaperState;
    mtPositive=supported&~conf&reviewOK&ratio>=m.evidenceRatio&wLow>0;
    mtPositive=minimumEvents(mtPositive,t,track,fsOut,m);
    multitaperState(mtPositive)=1; multitaperState(negative)=-1;
    oscPositive=valid&~conf&reviewOK&secondary.onScore>=m.oscillator.onThreshold;
    oscPositive=minimumEvents(oscPositive,t,secondary.frequency,fsOut,m);
    oscillatorState(oscPositive)=1; % No standalone oscillator amplitude-limit certificate.
    known=state~=0; frequency=NaN;
    if any(state==1), frequency=median(track(state==1)); end
    onPower=wPower; onPower(isfinite(onPower)&onPower<0)=0;
    active=state==1; ep=tbench.episodes(active,t,track,onPower,valid,0);
    for z=1:numel(ep), ep(z).boundaryUncertaintySeconds=m.windowSeconds/2+trial.clock.uncertaintySeconds; end
    row=struct('name',s.name,'kind',s.kind,'unit',s.unit,'frame',s.frame,'muscles',s.muscles, ...
        'time',t,'fs',fsOut,'coordinates',q,'slow',slow,'residual',y,'waveform',witness, ...
        'valid',valid,'state',state,'multitaperState',multitaperState,'oscillatorState',oscillatorState, ...
        'frequencyTrack',track,'frequency',frequency,'powerSigned',wPower, ...
        'powerLower',wLow,'powerUpper',wHigh,'meanBandPower',meanPower, ...
        'upperWindowPower',finiteQuantile(onPower,.95),'positiveSeconds',sum(state==1)/fsOut, ...
        'negativeSeconds',sum(state==-1)/fsOut,'indeterminateSeconds',sum(valid&state==0)/fsOut, ...
        'validFraction',min(1,sum(valid)/(fsOut*nominalSeconds)), ...
        'coordinateValidFraction',mean(valid),'captureFraction',min(1,numel(t)/(fsOut*nominalSeconds)), ...
        'recordedSeconds',numel(t)/fsOut,'validSeconds',sum(valid)/fsOut,'nominalSeconds',nominalSeconds, ...
        'stateCoverage',sum(known)/max(sum(valid),1), ...
        'occupancy',sum(active)/max(sum(known),1),'limitRMS',s.limitRMS, ...
        'negativeReferenceCalibrated',calibrated,'kineticReviewed',reviewOK, ...
        'qualitySupplied',s.qualitySupplied,'events',ep,'frames',frames,'secondary',secondary, ...
        'transfer',audit,'passHz',pass,'maximumPreservedHarmonic',floor(pass/m.searchBand(2)), ...
        'harmonicFundamentalHz',m.searchBand(2), ...
        'powerIntegrationBand',integrationBand, ...
        'referenceScope','model_supported_kinematic_rhythm_not_muscle_or_diagnosis_truth', ...
        'intervalScope','conditional_HAC_requires_whole_pipeline_calibration');
    if ~any(known), row.occupancy=NaN; end
    if isempty(targets), targets=row; else, targets(end+1)=row; end %#ok<AGROW>
end
ref=struct('schemaVersion',1,'status','ok','task',task,'clock',trial.clock,'targets',targets, ...
    'sourceSamplingHz',trial.fs,'missingSamples',trial.missingSamples,'policyHash',policy);
end
function [P,f]=vectorSpectrum(y,fs,seconds,overlap,NW,K,step)
N=round(seconds*fs); hop=max(1,round(N*(1-overlap))); nfft=2^nextpow2(max(N,fs/step));
f=(0:nfft/2)'*fs/nfft; P=zeros(size(f)); count=0; W=tbench.tapers(N,'dpss',NW,K);
scale=ones(size(f)); scale(2:end-1)=2;
for a=1:hop:size(y,1)-N+1
    v=y(a:a+N-1,:); if any(~isfinite(v(:))), continue; end
    v=bsxfun(@minus,v,mean(v,1)); sumPower=zeros(size(f));
    for d=1:size(v,2)
        z=fft(bsxfun(@times,v(:,d),W),nfft,1)/sqrt(fs);
        sumPower=sumPower+mean(abs(z(1:numel(f),:)).^2,2).*scale;
    end
    P=P+sumPower; count=count+1;
end
if count>0, P=P/count; else, P(:)=NaN; end
end
function good=minimumEvents(good,t,f,fs,m)
R=tbench.validRuns(good);
for k=1:size(R,1)
    ix=R(k,1):R(k,2); duration=numel(ix)/fs;
    if duration<m.minimumEpisodeSeconds||sum(f(ix))/fs<m.minimumCycles, good(ix)=false; end
end
end
function q=finiteQuantile(x,p)
x=x(isfinite(x)); q=NaN; if ~isempty(x), q=tbench.quantile(x,p); end
end
