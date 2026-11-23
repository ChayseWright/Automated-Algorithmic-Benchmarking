function rows = metrology(trial,mvc,pipelines,cfg)
%METROLOGY Complete-pipeline gain/polarity checks, plus reduction response.
ctx=tbench.preprocess.reference(trial,mvc,cfg); rows=struct([]);
for a=1:numel(pipelines)
    c=tbench.study.channel(ctx,1); base=execute(c,pipelines(a),cfg);
    for factor=[1.7 -1]
        changed=trial; changed.x=factor*trial.x;
        other=tbench.preprocess.reference(changed,abs(factor)*mvc,cfg);
        out=execute(tbench.study.channel(other,1),pipelines(a),cfg);
        ratio=NaN;
        if isfinite(base.powerMean)&&base.powerMean>0, ratio=out.powerMean/(factor^2*base.powerMean); end
        frequencyError=NaN;
        if isfinite(base.frequency)&&isfinite(out.frequency), frequencyError=abs(base.frequency-out.frequency); end
        row=struct('pipelineId',pipelines(a).id,'candidateId',pipelines(a).candidateId, ...
            'factor',factor,'baseStatus',base.status,'changedStatus',out.status,'powerGainRatio',ratio, ...
            'frequencyDifference',frequencyError,'detectionChanged',out.detected~=base.detected, ...
            'envelopePassErrorDB',ctx.reductionAudit.passErrorDB, ...
            'antiAliasStopAttenuationDB',ctx.reductionAudit.stopAttenuationDB, ...
            'originalFs',trial.fs,'analysisFs',ctx.fs,'carrierFs',ctx.reducedRawFs);
        if isempty(rows), rows=row; else, rows(end+1)=row; end %#ok<AGROW>
    end
end
end
function out=execute(c,p,cfg)
cfg=tbench.motion.inferenceConfig(cfg);
try
    pp=p.preprocess(c,cfg,p); out=p.analyze(pp,cfg,p);
catch err
    if cfg.failFast, rethrow(err); end
    out=tbench.emptyOutput(c); out.status='method_error'; out.message=err.message;
end
end
