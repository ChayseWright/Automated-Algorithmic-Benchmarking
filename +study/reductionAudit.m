function rows = reductionAudit(rates,cfg)
%REDUCTIONAUDIT G0 tones/harmonics, rejection and a pure-FM carrier challenge.
rows=struct([]); rates=rates(isfinite(rates)&rates>2*cfg.rawBand(2));
for fs=rates(:)'
    t=(0:ceil(12*fs)-1)'/fs;
    for f=[3 6 9 18 27 90]
        if f>=fs/2, continue; end
        x=cos(2*pi*f*t);
        [y,u,audit]=tbench.preprocess.rateChange(x,fs,cfg.analysisFs, ...
            cfg.reduction.envelopePassHz,cfg.reduction.envelopeStopHz,cfg.reduction);
        phi=2*pi*f*u; valid=isfinite(y); G=[ones(size(u)) cos(phi) sin(phi)];
        beta=pinv(G(valid,:))*y(valid); residual=y(valid)-G(valid,:)*beta;
        r=struct('kind','envelope_tone','inputFs',fs,'outputFs',cfg.analysisFs, ...
            'inputFrequencyHz',f,'modulationHz',NaN,'gain',norm(beta(2:3)), ...
            'phaseErrorRadians',atan2(-beta(3),beta(2)), ...
            'residualRMS',sqrt(mean(residual.^2)),'AAStopAttenuationDB',audit.stopAttenuationDB, ...
            'AAHalfSupportSeconds',audit.halfSupportSeconds,'validSamples',sum(valid));
        if isempty(rows), rows=r; else, rows(end+1)=r; end %#ok<AGROW>
    end
    target=fs; factor=max(1,floor(fs/cfg.reduction.minimumCarrierFs+1e-9));
    if cfg.reduction.enabled, target=fs/factor; end
    phase=2*pi*150*t+(10/6)*sin(2*pi*6*t); x=cos(phase);
    if abs(target-fs)<1e-9
        y=x; u=t; audit=struct('stopAttenuationDB',Inf,'halfSupportSeconds',0);
    else
        [y,u,audit]=tbench.preprocess.rateChange(x,fs,target,cfg.rawBand(2), ...
            min(target/2*.999,(cfg.rawBand(2)+target/2)/2),cfg.reduction);
    end
    expected=cos(2*pi*150*u+(10/6)*sin(2*pi*6*u)); good=isfinite(y);
    error=sqrt(mean((y(good)-expected(good)).^2));
    r.kind='raw_carrier_pure_FM'; r.outputFs=target; r.inputFrequencyHz=150; r.modulationHz=6;
    r.gain=NaN; r.phaseErrorRadians=NaN; r.residualRMS=error;
    r.AAStopAttenuationDB=audit.stopAttenuationDB; r.AAHalfSupportSeconds=audit.halfSupportSeconds;
    r.validSamples=sum(good); rows(end+1)=r; %#ok<AGROW>
end
end
