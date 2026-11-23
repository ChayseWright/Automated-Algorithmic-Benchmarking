function out = welchPower(c,cfg,frequency,detected,audit)
%WELCHPOWER Conventional signed peak-area estimate on the fixed reference PSD.
out=tbench.emptyOutput(c); out.audit=audit; out.valid=c.valid;
out.audit.proposalFrequency=frequency;
out.audit.powerEstimator='fixed_reference_Welch_signed_peak_area';
out.audit.reportingWindowSeconds=cfg.spectralWindow;
[P,f]=tbench.spectrum(c.envelope,c.fs,cfg.spectralWindow,cfg.spectralOverlap,'hann');
[background,~]=tbench.cfar(P,f,cfg);
if ~any(isfinite(P))||~isfinite(frequency), out.status='insufficient'; return; end
band=abs(f-frequency)<=cfg.peakHalfWidth&isfinite(background)&isfinite(P);
if sum(band)<2, out.status='insufficient'; return; end
out.signedExcessPower=trapz(f(band),P(band)-background(band));
out.powerMean=max(0,out.signedExcessPower); out.detected=logical(detected);
out.audit.reportingPSD=P; out.audit.reportingFrequency=f;
if out.detected, out.frequency=frequency; end
end
