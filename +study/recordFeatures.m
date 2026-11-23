function [row,events] = recordFeatures(out,c,cfg,meta)
%RECORDFEATURES Keep physical scalars and actual episode/window definitions.
[row,events]=tbench.summarize(out,c,cfg);
row.localPower50=NaN; row.localPower90=NaN; row.activeWindowPower90=NaN;
row.episodePowerMax=NaN; row.episodePower90=NaN; row.episodeDuration90=NaN;
row.episodeEnergy=NaN; row.trackChangeRate=NaN; row.harmonicRatio=NaN;
row.mvcActiveWindowPower90=NaN; row.mvcEpisodePowerMax=NaN;
row.mvcEnvelopeVariance=NaN; row.mvcLocalPower50=NaN; row.mvcLocalPower90=NaN;
row.dynamic=out.capabilities.dynamic; row.reviewFraction=0;
row.nativeDetectorPower=NaN; row.nativeDetectorFrequency=NaN; row.powerEstimator='unavailable';
row.tonicEnergy=NaN; row.mvcTonicEnergy=NaN;
if isfield(c,'energy')
    v=c.energy(c.valid&isfinite(c.energy));
    if ~isempty(v), row.tonicEnergy=mean(v); end
    if isfinite(c.mvc)&&c.mvc>0, row.mvcTonicEnergy=1e4*row.tonicEnergy/c.mvc^2; end
end
if isfield(out,'audit')
    row.nativeDetectorPower=tbench.option(out.audit,'nativeDetectorPower',out.powerMean);
    row.nativeDetectorFrequency=tbench.option(out.audit,'nativeDetectorFrequency',out.frequency);
    row.powerEstimator=tbench.option(out.audit,'powerEstimator','method_native');
end
if isfield(out,'audit')&&isfield(out.audit,'reviewMask'), row.reviewFraction=mean(out.audit.reviewMask); end
if strcmp(out.status,'ok')&&out.capabilities.waveform&&out.capabilities.dynamic
    valid=out.valid&isfinite(out.fundamental); active=logical(out.active)&valid;
    window=max(1,round(c.fs)); hop=max(1,round(c.fs/4)); values=[];
    R=tbench.validRuns(valid,window);
    for k=1:size(R,1)
        for start=R(k,1):hop:R(k,2)-window+1
            values(end+1)=mean(out.fundamental(start:start+window-1).^2); %#ok<AGROW>
        end
    end
    if ~isempty(values)
        row.localPower50=tbench.quantile(values,.5); row.localPower90=tbench.quantile(values,.9);
    end
    values=[]; R=tbench.validRuns(active);
    for k=1:size(R,1)
        if R(k,2)-R(k,1)+1<window
            values(end+1)=mean(out.fundamental(R(k,1):R(k,2)).^2); %#ok<AGROW>
        else
            for start=R(k,1):hop:R(k,2)-window+1
                values(end+1)=mean(out.fundamental(start:start+window-1).^2); %#ok<AGROW>
            end
        end
    end
    if ~isempty(values), row.activeWindowPower90=tbench.quantile(values,.9); end
    if ~isempty(events)
        row.episodePowerMax=max([events.power]); row.episodePower90=tbench.quantile([events.power],.9);
        row.episodeDuration90=tbench.quantile([events.duration],.9);
        row.episodeEnergy=sum([events.power].*[events.duration]);
    else, row.episodeEnergy=0;
    end
    if row.powerMean>0&&isfinite(row.harmonicPowerMean), row.harmonicRatio=row.harmonicPowerMean/row.powerMean; end
    ft=out.frequencyTrack; use=active(1:end-1)&active(2:end)&isfinite(ft(1:end-1))&isfinite(ft(2:end));
    delta=abs(diff(ft)); if any(use), row.trackChangeRate=mean(delta(use))*c.fs; end
end
if isfinite(c.mvc)&&c.mvc>0
    scale=1e4/c.mvc^2;
    row.mvcActiveWindowPower90=scale*row.activeWindowPower90;
    row.mvcEpisodePowerMax=scale*row.episodePowerMax;
    row.mvcEnvelopeVariance=scale*row.envelopeVariance;
    row.mvcLocalPower50=scale*row.localPower50; row.mvcLocalPower90=scale*row.localPower90;
end
names=fieldnames(meta); for k=1:numel(names), row.(names{k})=meta.(names{k}); end
for j=1:numel(events)
    for k=1:numel(names), events(j).(names{k})=meta.(names{k}); end
end
end
