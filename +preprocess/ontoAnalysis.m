function [z,audit] = ontoAnalysis(x,fs,c,cfg)
%ONTOANALYSIS Anti-alias a detector representation onto the complete clock.
p=c.preprocessAudit.parameters;
[v,t,audit]=tbench.preprocess.rateChange(x,fs,c.fs,p.envelopePassHz,p.envelopeStopHz,cfg.reduction);
index=round(t*c.fs)+1; inside=index>=1&index<=numel(c.t);
if any(abs(t(inside)-c.t(index(inside)))>1e-9), error('tbench:RepresentationClock','Detector and reporting clocks differ.'); end
z=nan(size(c.t)); z(index(inside))=v(inside);
end
