function out = panelSummary(M,cfg)
%PANELSUMMARY Subject-first mechanical panel with explicit evaluable counts.
out=struct('agreementSummary',struct([]),'waveformSummary',struct([]), ...
    'coverage',struct([]),'residualRepeatability',struct([]));
A=M.referenceAudit;
if ~isempty(A)
    for k=1:numel(A)
        r=A(k);
        out.coverage=append(out.coverage,struct('subject',r.subject,'group',r.group,'phenotype',r.phenotype, ...
            'task',r.task,'trial',r.trial,'target',r.target,'status',r.status, ...
            'validFraction',r.validFraction,'captureFraction',r.captureFraction, ...
            'recordedSeconds',r.recordedSeconds,'validSeconds',r.validSeconds,'nominalSeconds',r.nominalSeconds, ...
            'stateCoverage',r.stateCoverage, ...
            'positiveSeconds',r.positiveSeconds,'negativeSeconds',r.negativeSeconds, ...
            'indeterminateSeconds',r.indeterminateSeconds,'phaseReliable',r.phaseReliable, ...
            'timingReliable',r.timingReliable,'negativeReferenceCalibrated',r.negativeReferenceCalibrated, ...
            'referenceHash',r.referenceHash,'message',r.message));
    end
end
Q=M.agreement;
if ~isempty(Q)
    ids=unique({Q.pipelineId}); variants=unique({Q.referenceVariant}); targets=unique({Q.target});
    fields={'frequencyAbsError','frequencyMissing','trackingMissingFraction','frequencyWithinToleranceRecall', ...
        'recall','precisionResolved','IoU', ...
        'missedMotionBurden','falseAlarmSecondsPerMinute','onsetAbsError','offsetAbsError', ...
        'onsetBeyondUncertainty','algorithmCoverage','referenceValidFraction','stateCoverage'};
    for a=1:numel(ids)
        for b=1:numel(variants)
            for c=1:numel(targets)
                P=Q(strcmp({Q.pipelineId},ids{a})&strcmp({Q.referenceVariant},variants{b})&strcmp({Q.target},targets{c}));
                for subgroup={'all','control','ET','lower_severity','moderate_severity','higher_severity'}
                    R=subset(P,subgroup{1}); if isempty(R), continue; end
                    for type={'all','P','K'}
                        S=taskSubset(R,type{1}); if isempty(S), continue; end
                        r=struct('pipelineId',ids{a},'referenceVariant',variants{b},'target',targets{c}, ...
                            'subgroup',subgroup{1},'taskType',type{1},'subjects',numel(unique({S.subject})), ...
                            'trials',numel(S),'positiveTrials',sum([S.referencePositiveSeconds]>0), ...
                            'frequencyEligibleTrials',sum(isfinite([S.referenceFrequency])), ...
                            'negativeEligibleTrials',sum([S.referenceNegativeSeconds]>0), ...
                            'temporalCapableTrials',sum([S.nativeTemporalAvailable]));
                        for j=1:numel(fields)
                            key=fields{j}; T=S;
                            if any(strcmp(key,{'frequencyAbsError','frequencyMissing'})), T=T(isfinite([T.referenceFrequency])); end
                            [v,n]=macro(T,key,type{1},cfg); r.(key)=v; r.([key 'Subjects'])=n;
                        end
                        r.aggregation='trial_mean_then_task_mean_then_equal_PK_then_subject_mean';
                        out.agreementSummary=append(out.agreementSummary,r);
                    end
                end
            end
        end
    end
end
Q=M.waveformPrediction;
if ~isempty(Q)
    ids=unique({Q.pipelineId}); targets=unique({Q.target});
    for a=1:numel(ids)
        for b=1:numel(targets)
            P=Q(strcmp({Q.pipelineId},ids{a})&strcmp({Q.target},targets{b}));
            for subgroup={'all','control','ET','lower_severity','moderate_severity','higher_severity'}
                R=subset(P,subgroup{1}); if isempty(R), continue; end
                r=struct('pipelineId',ids{a},'target',targets{b},'subgroup',subgroup{1}, ...
                    'subjects',numel(unique({R.subject})),'trials',numel(R),'evaluableTrials',sum(strcmp({R.status},'ok')), ...
                    'scope','person_calibrated_repetition_validation');
                fields={'gainOverNull','gainOverReference','inputCoverage','frozenTemporalMeanRetained', ...
                    'F0waveformRetained','F1waveformRetained','F2waveformRetained','F3waveformRetained','RsubjectWaveformRetained'};
                for j=1:numel(fields)
                    key=fields{j}; [v,n]=macro(R,key,'all',cfg); r.(key)=v; r.([key 'Subjects'])=n;
                end
                out.waveformSummary=append(out.waveformSummary,r);
            end
        end
    end
end
Q=M.predictions;
if ~isempty(Q)
    values=arrayfun(@(r)[r.pipelineId '|' r.level '|' r.target '|' r.subject '|' r.task],Q,'UniformOutput',false);
    keys=unique(values);
    for k=1:numel(keys)
        R=Q(strcmp(values,keys{k})); if numel(R)<2, continue; end
        e=[R.predictionTransformed]-[R.referenceTransformed]; sd=std(e,0);
        % Repeatability of error conditional on recorded motion, not raw EMG ICC.
        r=struct('pipelineId',R(1).pipelineId,'level',R(1).level,'target',R(1).target,'subject',R(1).subject, ...
            'group',R(1).group,'severityGroup',R(1).severityGroup,'task',R(1).task,'repetitions',numel(R), ...
            'meanTransformedError',mean(e),'withinTaskErrorSD',sd,'scope','motion_conditioned_decoder_residual');
        out.residualRepeatability=append(out.residualRepeatability,r);
    end
end
end
function [value,n]=macro(rows,key,type,cfg)
value=NaN; n=0; if isempty(rows), return; end
ids=unique({rows.subject}); allValues=[];
for s=1:numel(ids)
    Q=rows(strcmp({rows.subject},ids{s})); kinds=type;
    if strcmp(type,'all'), kinds='PK'; end
    v=[];
    for kind=kinds
        tasks=cfg.tasks(cellfun(@(x)x(1)==kind,cfg.tasks)); taskValues=[];
        for j=1:numel(tasks)
            R=Q(strcmp({Q.task},tasks{j})); if isempty(R), continue; end
            z=[R.(key)]; z=z(isfinite(z)); if ~isempty(z), taskValues(end+1)=mean(z); end %#ok<AGROW>
        end
        if ~isempty(taskValues), v(end+1)=mean(taskValues); end %#ok<AGROW>
    end
    if numel(v)==numel(kinds), allValues(end+1)=mean(v); end %#ok<AGROW>
end
n=numel(allValues); if n>0, value=mean(allValues); end
end
function r=taskSubset(r,type)
if ~strcmp(type,'all'), r=r(cellfun(@(x)x(1)==type,{r.task})); end
end
function rows=subset(rows,name)
if isempty(rows), return; end
switch name
    case 'control', rows=rows(strcmp({rows.group},'control'));
    case 'ET', rows=rows(strcmp({rows.group},'ET'));
    case 'lower_severity', rows=rows(strcmp({rows.group},'ET')&[rows.severityGroup]==1);
    case 'moderate_severity', rows=rows(strcmp({rows.group},'ET')&[rows.severityGroup]==2);
    case 'higher_severity', rows=rows(strcmp({rows.group},'ET')&[rows.severityGroup]==3);
end
end
function a=append(a,b)
if isempty(a), a=b(:); else, a=[a(:);b(:)]; end
end
