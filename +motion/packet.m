function p = packet(p,out,m,c)
%PACKET Temporary one-trial phase-preserving output; never retained per candidate.
if isempty(p)
    N=numel(c.t);
    p=struct('time',c.t,'fs',c.fs,'waveform',nan(N,15),'valid',false(N,15), ...
        'active',false(N,15),'frequencyTrack',nan(N,15),'frequency',nan(1,15), ...
        'dynamic',false(1,15),'hasWaveform',false(1,15),'records',struct([]));
end
if numel(out.valid)~=numel(c.t), return; end
p.valid(:,m)=logical(out.valid)&strcmp(out.status,'ok'); p.frequency(m)=out.frequency;
p.dynamic(m)=out.capabilities.dynamic; p.hasWaveform(m)=out.capabilities.waveform;
if out.capabilities.dynamic&&~isempty(out.active)
    p.active(:,m)=logical(out.active)&p.valid(:,m); p.frequencyTrack(:,m)=out.frequencyTrack;
end
if out.capabilities.waveform&&~isempty(out.fundamental), p.waveform(:,m)=out.fundamental; end
end
