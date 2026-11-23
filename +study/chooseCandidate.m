function [index,status] = chooseCandidate(scores,cfg)
%CHOOSECANDIDATE Mild capture first, subject to control, severe and fidelity gates.
if isempty(scores), error('tbench:Selection','No candidates were evaluated.'); end
eligible=[scores.eligible]; status='constraints_met';
if ~any(eligible)
    % Retain the measured failure instead of silently relaxing the constraints.
    status='constraints_unmet'; violation=[scores.violation];
    violation(~isfinite(violation))=realmax; best=min(violation);
    eligible=abs(violation-best)<=cfg.tuning.tieTolerance;
end
ix=find(eligible); objective=nan(size(ix));
for k=1:numel(ix)
    m=scores(ix(k)).metrics;
    if strcmp(cfg.clinical.priority,'severe'), objective(k)=m.severeCapture; else, objective(k)=m.mildCapture; end
end
objective(~isfinite(objective))=-Inf; maximum=max(objective); ix=ix(objective>=maximum-cfg.tuning.tieTolerance);
if isempty(ix), ix=find(eligible); end
if isfield(cfg,'motion')&&cfg.motion.enabled&&cfg.motion.selectionTieBreak&&isfield(scores,'motion')
    gain=arrayfun(@(s)s.motion.gain,scores(ix));
    if all(isfinite(gain))
        best=max(gain); ix=ix(gain>=best-cfg.tuning.tieTolerance);
    end
end
loss=arrayfun(@(s)s.metrics.balancedLogLoss,scores(ix)); loss(~isfinite(loss))=Inf;
minimum=min(loss); ix=ix(loss<=minimum+cfg.tuning.tieTolerance);
if isempty(ix), ix=find(eligible); end
% Deterministic final tie-breaker favors the prespecified simpler first entry.
index=ix(1);
end
