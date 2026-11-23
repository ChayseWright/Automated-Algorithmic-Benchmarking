function [summary,contrasts,retention] = clinicalSummary(P,cfg)
%CLINICALSUMMARY Paired participant inference from untouched outer predictions.
summary=struct([]); contrasts=struct([]); retention=struct([]);
if isempty(P), return; end
ids=unique({P.pipelineId}); levels=unique({P.level});
for a=1:numel(ids)
    for l=1:numel(levels)
        Q=P(strcmp({P.pipelineId},ids{a})&strcmp({P.level},levels{l}));
        if isempty(Q), continue; end
        m=tbench.benchmark.metrics(Q); m.pipelineId=ids{a}; m.level=levels{l};
        summary=append(summary,m);
    end
end
% Contrasts use a common participant set and are conditional on fitted folds.
for a=1:numel(ids)-1
    for b=a+1:numel(ids)
        A=P(strcmp({P.pipelineId},ids{a})&strcmp({P.level},cfg.clinical.featureLevel));
        B=P(strcmp({P.pipelineId},ids{b})&strcmp({P.level},cfg.clinical.featureLevel));
        [A,B]=align(A,B); if isempty(A), continue; end
        state=rng; rng(cfg.seed+901,'twister');
        d=zeros(cfg.tuning.bootstrapSamples,6); strata={find(strcmp({A.group},'control'))};
        for g=0:3, strata{end+1}=find(strcmp({A.group},'ET')&[A.trueSeverityGroup]==g); end %#ok<AGROW>
        for z=1:size(d,1)
            ix=[];
            for g=1:numel(strata)
                v=strata{g}; if ~isempty(v), ix=[ix v(randi(numel(v),1,numel(v)))]; end %#ok<AGROW>
            end
            ma=tbench.benchmark.metrics(A(ix)); mb=tbench.benchmark.metrics(B(ix));
            d(z,:)=[ma.mildCapture-mb.mildCapture ma.controlFPR-mb.controlFPR ma.balancedLogLoss-mb.balancedLogLoss ...
                ma.severeCapture-mb.severeCapture ma.falseSevereLower-mb.falseSevereLower ma.clinicalCoverage-mb.clinicalCoverage];
        end
        rng(state); ma=tbench.benchmark.metrics(A); mb=tbench.benchmark.metrics(B);
        q1=tbench.quantile(d(:,1),[.025 .975]); q2=tbench.quantile(d(:,2),[.025 .975]); q3=tbench.quantile(d(:,3),[.025 .975]);
        q4=tbench.quantile(d(:,4),[.025 .975]); q5=tbench.quantile(d(:,5),[.025 .975]); q6=tbench.quantile(d(:,6),[.025 .975]);
        r=struct('pipelineA',ids{a},'pipelineB',ids{b},'participants',numel(A), ...
            'mildCaptureDifference',ma.mildCapture-mb.mildCapture,'mildCI025',q1(1),'mildCI975',q1(2), ...
            'controlFPRDifference',ma.controlFPR-mb.controlFPR,'controlCI025',q2(1),'controlCI975',q2(2), ...
            'logLossDifference',ma.balancedLogLoss-mb.balancedLogLoss,'lossCI025',q3(1),'lossCI975',q3(2), ...
            'severeCaptureDifference',ma.severeCapture-mb.severeCapture,'severeCI025',q4(1),'severeCI975',q4(2), ...
            'falseSevereLowerDifference',ma.falseSevereLower-mb.falseSevereLower,'falseSevereCI025',q5(1),'falseSevereCI975',q5(2), ...
            'coverageDifference',ma.clinicalCoverage-mb.clinicalCoverage,'coverageCI025',q6(1),'coverageCI975',q6(2), ...
            'uncertaintyScope','paired_subject_bootstrap_conditional_on_outer_fits');
        contrasts=append(contrasts,r);
    end
end
for a=1:numel(ids)
    rich=P(strcmp({P.pipelineId},ids{a})&strcmp({P.level},'R'));
    if isempty(rich), continue; end
    for l=1:numel(levels)
        if strcmp(levels{l},'R'), continue; end
        Q=P(strcmp({P.pipelineId},ids{a})&strcmp({P.level},levels{l}));
        [Q,R]=align(Q,rich); if isempty(Q), continue; end
        m=tbench.benchmark.metrics(Q); mr=tbench.benchmark.metrics(R);
        base=log(2); gain=base-mr.balancedLogLoss; retained=NaN;
        if isfinite(gain)&&gain>1e-8, retained=(base-m.balancedLogLoss)/gain; end
        r=struct('pipelineId',ids{a},'level',levels{l},'participants',numel(Q), ...
            'baselineBalancedLogLoss',base,'richBalancedLogLoss',mr.balancedLogLoss, ...
            'compressedBalancedLogLoss',m.balancedLogLoss, ...
            'logLossIncrease',m.balancedLogLoss-mr.balancedLogLoss, ...
            'retainedDiagnosticUtility',retained,'mildCaptureChange',m.mildCapture-mr.mildCapture, ...
            'controlFPRChange',m.controlFPR-mr.controlFPR,'scope','operational_decoder_not_total_information');
        if strcmp(levels{l},'B0'), r.scope='alternative_tonic_MVC_baseline_not_nested_compression'; end
        retention=append(retention,r);
    end
end
end
function [A,B]=align(A,B)
if isempty(A)||isempty(B), A=[]; B=[]; return; end
[~,ia,ib]=intersect({A.subject},{B.subject},'stable'); A=A(ia); B=B(ib);
end
function a=append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
