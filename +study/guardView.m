function row = guardView(c,cfg)
%GUARDVIEW Fixed broad-band witness for rhythms outside the primary search.
row=struct('guardPeakHz',NaN,'guardStatistic',NaN,'guardAccepted',false,'guardOutsidePrimary',false);
[P,f]=tbench.spectrum(c.envelope,c.fs,cfg.spectralWindow,cfg.spectralOverlap,'hann');
if ~any(isfinite(P)), return; end
[~,ratio]=tbench.cfar(P,f,cfg);
ix=find(f>=cfg.guardBand(1)&f<=cfg.guardBand(2)&isfinite(ratio));
if isempty(ix), return; end
[value,j]=max(ratio(ix)); row.guardPeakHz=f(ix(j)); row.guardStatistic=value;
row.guardAccepted=value>=cfg.cfar.ratio;
row.guardOutsidePrimary=row.guardAccepted&&(row.guardPeakHz<cfg.searchBand(1)||row.guardPeakHz>cfg.searchBand(2));
end
