function gate = calibrationGate(rows,cfg)
%CALIBRATIONGATE Development challenge thresholds, with explicit case counts.
gate=struct('passed',true,'coverage',NaN,'strongCases',0,'strongDetection',NaN, ...
    'strongPowerError',NaN,'strongFrequencyError',NaN,'nullCases',0,'nullAlarmFraction',NaN, ...
    'weakCases',0,'weakFrequencyRecovered',NaN,'lateStrongCases',0,'lateStrongRecovered',NaN,'violation',0);
if isempty(rows)
    gate.passed=~cfg.tuning.calibrationGate; if ~gate.passed, gate.violation=1; end
    return;
end
ok=strcmp({rows.status},'ok'); gate.coverage=mean(ok);
strong=[rows.depth]>=.6&[rows.intermittent]; gate.strongCases=sum(strong);
gate.nullCases=sum([rows.depth]==0);
if gate.nullCases>0, gate.nullAlarmFraction=mean([rows([rows.depth]==0).detected]); end
if any(strong)
    gate.strongDetection=sum([rows(strong).detected]&ok(strong))/sum(strong);
    errors=[rows(strong).relativePowerError]; errors(~isfinite(errors))=Inf;
    gate.strongPowerError=median(errors);
    errors=[rows(strong).frequencyAbsError]; errors(~isfinite(errors))=Inf;
    gate.strongFrequencyError=median(errors);
end
weak=[rows.depth]>0&[rows.depth]<=.10; gate.weakCases=sum(weak);
if any(weak), gate.weakFrequencyRecovered=mean([rows(weak).detected]&ok(weak)&[rows(weak).frequencyAbsError]<=.75); end
if isfield(rows,'pattern')
    late=strong&strcmp({rows.pattern},'late'); gate.lateStrongCases=sum(late);
    if any(late), gate.lateStrongRecovered=mean([rows(late).detected]&ok(late)&[rows(late).frequencyAbsError]<=.75); end
end
if cfg.tuning.calibrationGate
    if ~isfinite(gate.coverage), gate.violation=gate.violation+1;
    else, gate.violation=gate.violation+max(0,cfg.tuning.minCalibrationCoverage-gate.coverage);
    end
    if gate.strongCases==0
        gate.violation=gate.violation+1;
    else
        gate.violation=gate.violation+max(0,cfg.tuning.minStrongDetection-gate.strongDetection);
        if ~isfinite(gate.strongPowerError), gate.violation=gate.violation+1;
        else, gate.violation=gate.violation+max(0,gate.strongPowerError-cfg.tuning.maxStrongPowerError); end
        if ~isfinite(gate.strongFrequencyError), gate.violation=gate.violation+1;
        else, gate.violation=gate.violation+max(0,gate.strongFrequencyError-cfg.tuning.maxStrongFrequencyError); end
    end
    if gate.nullCases==0, gate.violation=gate.violation+1;
    else, gate.violation=gate.violation+max(0,gate.nullAlarmFraction-cfg.tuning.maxNullAlarmFraction); end
    gate.passed=gate.violation<=1e-12;
end
end
