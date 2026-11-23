function prepared = detector(c,cfg,p)
%DETECTOR Shared detector-domain filtering; reporting reference stays fixed.
prepared=c; prepared.referenceEnvelope=c.envelope; prepared.referenceY=c.y;
prepared.referenceEnergy=c.energy; prepared.referenceValid=c.valid;
prepared.referenceDynamicValid=c.dynamicValid; prepared.detectorDomainChanged=false;
if p.rawHighpassHz~=cfg.rawBand(1)||p.envelopePassHz~=cfg.reduction.envelopePassHz||p.envelopeStopHz~=cfg.reduction.envelopeStopHz
    band=[p.rawHighpassHz cfg.rawBand(2)];
    carrier=tbench.applyFIR(c.raw,tbench.designFIR(c.rawFs,band,cfg.carrierFIRSeconds));
    [e,~,~]=tbench.preprocess.rateChange(abs(carrier),c.rawFs,c.fs,p.envelopePassHz,p.envelopeStopHz,cfg.reduction);
    [q,~,~]=tbench.preprocess.rateChange(carrier.^2,c.rawFs,c.fs,p.envelopePassHz,p.envelopeStopHz,cfg.reduction);
    q(q<=0)=NaN; valid=isfinite(e)&isfinite(q)&c.valid; e(~valid)=NaN; q(~valid)=NaN;
    slow=tbench.applyFIR(e,tbench.designFIR(c.fs,cfg.baselineCutoff,cfg.baselineFIRSeconds));
    prepared.envelope=e; prepared.energy=q; prepared.y=e-slow;
    prepared.valid=valid; prepared.dynamicValid=valid&isfinite(prepared.y); prepared.detectorDomainChanged=true;
end
end
