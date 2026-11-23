function r = scoreTrial(packet,ref,target,m,declaredTemporal)
%SCORETRIAL Mechanical association; absence of motion is not a muscle label.
if nargin<5, declaredTemporal=~isempty(packet)&&all(packet.dynamic(target.muscles)); end
r=struct('target',target.name,'kind',target.kind,'unit',target.unit,'status','ok', ...
    'referenceVariant','consensus', ...
    'frequencyScope','trial_level_reference_supported','frequencyAbsError',NaN,'frequencyMissing',true, ...
    'acceptedTrackSeconds',NaN,'trackingMissingFraction',NaN,'frequencyWithinToleranceRecall',NaN, ...
    'harmonicRelated',false,'referenceFrequency',target.frequency,'estimatedFrequency',NaN, ...
    'referencePositiveSeconds',target.positiveSeconds,'referenceNegativeSeconds',target.negativeSeconds, ...
    'referenceIndeterminateSeconds',target.indeterminateSeconds,'referenceValidFraction',target.validFraction, ...
    'stateCoverage',target.stateCoverage,'nativeTemporalAvailable',false,'algorithmCoverage',NaN, ...
    'TPseconds',NaN,'FPseconds',NaN,'FNseconds',NaN,'recall',NaN,'precisionResolved',NaN,'IoU',NaN, ...
    'unconfirmedPositiveSeconds',NaN,'missedMotionBurden',NaN,'falseAlarmSecondsPerMinute',NaN, ...
    'referenceEpisodes',numel(target.events),'matchedEpisodes',NaN,'onsetAbsError',NaN,'offsetAbsError',NaN, ...
    'onsetBeyondUncertainty',NaN,'offsetBeyondUncertainty',NaN, ...
    'boundaryUncertaintySeconds',m.windowSeconds/2+ref.clock.uncertaintySeconds, ...
    'negativeReferenceCalibrated',target.negativeReferenceCalibrated,'timingReliable',ref.clock.timingReliable);
r.nativeTemporalAvailable=declaredTemporal;
sourceClockOK=numel(target.muscles)==1||tbench.option(m,'sourceSynchronized',false);
if ~sourceClockOK, r.timingReliable=false; end
if isempty(packet)
    r.status='emg_unavailable';
    if declaredTemporal&&ref.clock.timingReliable&&sourceClockOK
        r.algorithmCoverage=0;
        r.TPseconds=0; r.FPseconds=0; r.FNseconds=target.positiveSeconds;
        if target.positiveSeconds>0
            r.recall=0; r.IoU=0; r.missedMotionBurden=1; r.acceptedTrackSeconds=0;
            r.trackingMissingFraction=1; r.frequencyWithinToleranceRecall=0;
        end
        if target.negativeSeconds>0, r.falseAlarmSecondsPerMinute=0; end
        r.unconfirmedPositiveSeconds=0; r.matchedEpisodes=0;
    end
    return;
end
M=target.muscles; f=packet.frequency(M); f=f(isfinite(f));
if ~isempty(f), r.estimatedFrequency=median(f); r.frequencyMissing=false; end
if isfinite(target.frequency)&&isfinite(r.estimatedFrequency)
    r.frequencyAbsError=abs(target.frequency-r.estimatedFrequency);
    ratio=r.estimatedFrequency/target.frequency;
    r.harmonicRelated=min(abs(ratio-[.5 2]))<.10;
end
if ~r.nativeTemporalAvailable, r.status='native_temporal_unavailable'; return; end
if ~sourceClockOK, r.status='EMG_source_clock_unverified'; return; end
if ~ref.clock.timingReliable, r.status='timing_unverified_or_uncertain'; return; end
nearest=interp1(packet.time,(1:numel(packet.time))',target.time,'nearest',NaN); inside=isfinite(nearest);
active=false(size(target.time)); supported=false(size(active));
ix=nearest(inside); V=packet.valid(ix,M);
supported(inside)=mean(V,2)>=m.minimumSourceCoverage;
active(inside)=sum(packet.active(ix,M)&V,2)>=m.minimumActiveMuscles;
active=active&supported;
valid=target.valid; positive=target.state==1&valid; negative=target.state==-1&valid; unknown=target.state==0&valid;
r.algorithmCoverage=sum(supported&valid)/max(sum(valid),1);
tp=sum(active&positive); fn=sum(~active&positive); fp=sum(active&negative);
r.TPseconds=tp/target.fs; r.FPseconds=fp/target.fs; r.FNseconds=fn/target.fs;
r.recall=divide(tp,tp+fn); r.precisionResolved=divide(tp,tp+fp); r.IoU=divide(tp,tp+fn+fp);
r.unconfirmedPositiveSeconds=sum(active&unknown)/target.fs;
p=target.powerSigned; p(isfinite(p)&p<0)=0; good=positive&isfinite(p);
r.missedMotionBurden=divide(sum(p(good&~active)),sum(p(good)));
r.falseAlarmSecondsPerMinute=60*divide(fp,sum(negative));
if any(positive)
    track=nan(size(active)); track(positive)=target.frequencyTrack(positive);
    fnew=packet.frequencyTrack; aligned=nan(size(active));
    for k=find(inside(:))'
        v=fnew(nearest(k),M); v=v(isfinite(v)); if ~isempty(v), aligned(k)=median(v); end
    end
    both=positive&active&isfinite(aligned)&isfinite(track);
    r.acceptedTrackSeconds=sum(both)/target.fs;
    r.trackingMissingFraction=1-sum(both)/sum(positive);
    r.frequencyWithinToleranceRecall=sum(both&abs(aligned-track)<=m.frequencyAgreementHz)/sum(positive);
    if any(both), r.frequencyAbsError=median(abs(aligned(both)-track(both))); r.frequencyScope='detected_overlap_conditional'; end
end
found=tbench.validRuns(active&valid); truth=tbench.validRuns(positive); matches=[];
for a=1:size(truth,1)
    for b=1:size(found,1)
        overlap=max(0,min(truth(a,2),found(b,2))-max(truth(a,1),found(b,1))+1);
        union=max(truth(a,2),found(b,2))-min(truth(a,1),found(b,1))+1;
        if overlap>0, matches(end+1,:)=[overlap/union a b]; end %#ok<AGROW>
    end
end
usedTruth=[]; usedFound=[]; onset=[]; offset=[];
if ~isempty(matches)
    [~,order]=sort(matches(:,1),'descend');
    for index=order(:)'
        a=matches(index,2); b=matches(index,3);
        if any(usedTruth==a)||any(usedFound==b), continue; end
        usedTruth(end+1)=a; usedFound(end+1)=b; %#ok<AGROW>
        onset(end+1)=abs(found(b,1)-truth(a,1))/target.fs; offset(end+1)=abs(found(b,2)-truth(a,2))/target.fs; %#ok<AGROW>
    end
end
r.matchedEpisodes=numel(usedTruth);
if ~isempty(onset)
    r.onsetAbsError=median(onset); r.offsetAbsError=median(offset);
    r.onsetBeyondUncertainty=max(0,r.onsetAbsError-r.boundaryUncertaintySeconds);
    r.offsetBeyondUncertainty=max(0,r.offsetAbsError-r.boundaryUncertaintySeconds);
end
end
function q=divide(a,b)
q=NaN; if b>0, q=a/b; end
end
