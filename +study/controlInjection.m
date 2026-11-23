function rows = controlInjection(trial,mvc,pipelines,cfg,outputs,ctx)
%CONTROLINJECTION Held-out control counterfactual in the fixed measured domain.
% Truth is the known-phase projection of the EXACT reference-envelope change.
% It is not latent neural-drive power or an additive-power assumption.
rows=struct([]);
algorithmCfg=tbench.motion.inferenceConfig(cfg);
for m=cfg.injection.muscles(:)'
    c0=tbench.study.channel(ctx,m);
    for f=cfg.injection.frequencies(:)'
        for depth=cfg.injection.depths(:)'
            start=cfg.injection.startFraction*trial.duration;
            stop=min(trial.duration,start+cfg.injection.durationSeconds);
            gate=double(trial.t>=start&trial.t<stop);
            changed=trial; changed.x=trial.x(:,m).*(1+depth*gate.*cos(2*pi*f*trial.t));
            reference=tbench.preprocess.reference(changed,mvc(m),cfg); c1=tbench.study.channel(reference,1);
            valid=c0.dynamicValid&c1.dynamicValid;
            g=double(c0.t>=start&c0.t<stop); phi=2*pi*f*c0.t;
            B=[ones(size(phi)) g g.*cos(phi) g.*sin(phi)];
            delta=c1.y-c0.y; truth=nan(size(delta)); truthPower=NaN;
            if sum(valid)>=20
                beta=pinv(B(valid,:))*delta(valid); truth=B(:,3:4)*beta(3:4);
                truthPower=mean(truth(valid).^2);
            end
            for a=1:numel(pipelines)
                comparisonTruthPower=truthPower;
                try
                    pp=pipelines(a).preprocess(c1,algorithmCfg,pipelines(a)); out=pipelines(a).analyze(pp,algorithmCfg,pipelines(a));
                catch err
                    if cfg.failFast, rethrow(err); end
                    out=tbench.emptyOutput(c1); out.status='method_error'; out.message=err.message;
                end
                base=outputs{a,m}; incremental=NaN; relative=NaN; waveformAvailable=false;
                if strcmp(base.status,'ok')&&strcmp(out.status,'ok')&&base.capabilities.waveform&&out.capabilities.waveform
                    ok=valid&base.valid&out.valid&isfinite(truth);
                    if sum(ok)>=20
                        difference=out.fundamental-base.fundamental;
                        b=pinv(B(ok,:))*difference(ok); projected=B(:,3:4)*b(3:4);
                        incremental=mean(projected(ok).^2); comparisonTruthPower=mean(truth(ok).^2);
                        waveformAvailable=true;
                        if comparisonTruthPower>0, relative=abs(incremental-comparisonTruthPower)/comparisonTruthPower; end
                    end
                end
                recovered=out.detected&&isfinite(out.frequency)&&abs(out.frequency-f)<=.75;
                r=struct('pipelineId',pipelines(a).id,'candidateId',pipelines(a).candidateId, ...
                    'muscle',m,'frequency',f,'depth',depth,'status',out.status, ...
                    'truthScope','counterfactual_reference_change_known_phase_projection', ...
                    'truthPower',comparisonTruthPower,'incrementalProjectedPower',incremental, ...
                    'relativePowerError',relative,'waveformComparisonAvailable',waveformAvailable, ...
                    'frequencyRecovered',recovered,'estimatedFrequency',out.frequency, ...
                    'baseDetected',base.detected,'nativeTotalPowerChange',out.powerMean-base.powerMean);
                if isempty(rows), rows=r; else, rows(end+1)=r; end %#ok<AGROW>
            end
        end
    end
end
end
