function a = trialFeatures(F,muscles)
%TRIALFEATURES Actual per-record descriptor compression; no target information.
F=F(ismember([F.muscle],muscles)); good=strcmp({F.status},'ok')&isfinite([F.mvcPowerMean]);
a=struct('F0',nan(1,2),'F1',nan(1,3),'F2',nan(1,6),'F3',nan(1,9),'R',nan(1,15), ...
    'B0',nan(1,2),'coverage',0);
if isempty(F), return; end
a.coverage=sum([F(good).validFraction])/numel(muscles);
v=[F.mvcTonicEnergy]; mvc=[F.mvcReference];
mvc(~isfinite(mvc)|mvc<=0)=NaN; a.B0=[finite(v) finite(mvc)];
if ~any(good), return; end
G=F(good); power=[G.mvcPowerMean]; variance=finite([G.mvcEnvelopeVariance]);
base=[mean(power) max(power) mean([G.detected])]; activity=[G.occupancy]>0&isfinite([G.occupancy]);
occupancy=finite([G.occupancy]); if ~any([G.dynamic]), occupancy=0; end
extra=[occupancy off([G.mvcActiveWindowPower90],activity) mean(activity) ...
    off([G.mvcEpisodePowerMax],activity) off([G.episodeDuration90],activity) ...
    sumFinite([G.episodeCount])/max(sum([G.validSeconds])/60,eps) ...
    off([G.mvcLocalPower50],true(size(activity))) off([G.mvcLocalPower90],true(size(activity))) ...
    off([G.harmonicRatio],activity) off([G.frequencyIQR],activity) off([G.trackChangeRate],activity)];
if ~any([G.dynamic]), extra=zeros(1,11); end
a.F0=[variance base(3)]; a.F1=base; a.F2=[base extra(1:3)];
a.F3=[base extra(1:6)]; a.R=[base extra variance];
end
function v=finite(x)
x=x(isfinite(x)); v=NaN; if ~isempty(x), v=mean(x); end
end
function v=off(x,active)
if ~any(active), v=0; return; end
x=x(active&isfinite(x)); v=NaN; if ~isempty(x), v=tbench.quantile(x,.9); end
end
function v=sumFinite(x)
x=x(isfinite(x)); v=0; if ~isempty(x), v=sum(x); end
end
