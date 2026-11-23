function rows = calibrate(m,seed)
%CALIBRATE Independent synthetic rigid translations, not device calibration.
state=rng; cleanup=onCleanup(@()rng(state)); rng(seed+6817,'twister'); %#ok<NASGU>
rows=struct([]); fs=m.calibration.fs; T=m.calibration.duration; t=(0:round(fs*T)-1)'/fs;
caseIndex=0;
for pattern={'continuous','late'}
    for frequency=m.calibration.frequencies
        for amplitude=m.calibration.amplitudesRMS
            caseIndex=caseIndex+1; gate=ones(size(t));
            if strcmp(pattern{1},'late'), gate=double(t>=.70*T&t<.90*T); end
            direction=[1 .3 -.2]; direction=direction/norm(direction);
            signal=sqrt(2)*amplitude*bsxfun(@times,gate.*cos(2*pi*frequency*t),direction);
            slow=[.2*sin(2*pi*.25*t) zeros(numel(t),2)];
            values=signal+slow+m.calibration.noiseSD*randn(numel(t),3); valid=true(size(values));
            gap=mod(caseIndex,7)==0;
            if gap, valid(t>.4*T&t<.4*T+.15,:)=false; end
            target=struct('name','synthetic_translation','kind','translation','unit','mm','frame','calibrated_bench', ...
                'values',values,'valid',valid,'noiseSD',m.calibration.noiseSD,'limitRMS',3*m.calibration.noiseSD, ...
                'calibrationVerified',true,'calibrationPolicyHash',tbench.motion.policyHash(m),'reviewed',true);
            record=struct('time',t,'targets',target,'sync',struct('verified',true,'offsetSeconds',0,'rateRatio',1,'uncertaintySeconds',0));
            ref=tbench.motion.reference(record,struct('timeOrigin',0),m,'P1'); a=ref.targets(1);
            actual=interp1(t,double(gate),a.time,'nearest',0)>.5&amplitude>0;
            use=a.valid; accepted=a.state==1; recall=NaN;
            if any(actual&use), recall=sum(accepted&actual&use)/sum(actual&use); end
            ferror=NaN; if amplitude>0&&isfinite(a.frequency), ferror=abs(a.frequency-frequency); end
            row=struct('caseIndex',caseIndex,'scope','synthetic_motion_not_instrument_calibration', ...
                'pattern',pattern{1},'frequency',frequency,'truthActiveRMS',amplitude,'gap',gap, ...
                'referenceFrequency',a.frequency,'frequencyAbsError',ferror,'activeTimeRecall',recall, ...
                'falsePositiveSeconds',sum(accepted&~actual&use)/a.fs,'stateCoverage',a.stateCoverage, ...
                'validFraction',a.validFraction,'estimatedUpperPower',a.upperWindowPower, ...
                'measuredBandPower',a.meanBandPower,'referencePolicyHash',ref.policyHash);
            if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
        end
    end
end
end
