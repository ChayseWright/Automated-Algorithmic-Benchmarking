function [pp,ap,defaults,description] = parameterSpace(p,cfg)
%PARAMETERSPACE Method-specific parameters; benchmark definitions stay frozen.
C=tbench.methodCatalog(); j=find(strcmp({C.id},p.analysisId),1);
if isempty(j), error('tbench:Space','Unknown analysis ID.'); end
family=C(j).family; kind=C(j).variant;
pp=struct('rawHighpassHz',[20 35 40]); ap=struct(); defaults=struct();
description=struct('analysisId',p.analysisId,'family',family,'variant',kind, ...
    'strategy','prespecified_parameter_coverage_then_training_only_coordinate_refinement', ...
    'optimality','best_evaluated_configuration_in_declared_discrete_space', ...
    'reportingDomain','fixed_original_rate_reference_observer', ...
    'mildPriority',strcmp(cfg.clinical.priority,'mild'));
threshold=true; background=false; local=false;
switch family
    case 'legacy'
        if strcmp(kind,'welch')
            pp.windowSeconds=[2 3 4 6]; pp.spectralOverlap=[0 .5 .75]; background=true;
        else
            local=true; pp.localWindowSeconds=[1 1.5 2 3 4]; pp.localHop=[.1 .2 .3];
            ap.minCycles=[3 4 5]; ap.switchPenalty=[.5 1 2]; ap.frequencyPenalty=[.15 .35 .7]; background=true;
            if any(strcmp(kind,{'dpss','kalman'}))
                ap.NW=[1.5 2 2.5 3]; ap.nTapers=[1 2 3 4];
                defaults.NW=cfg.local.NW; defaults.nTapers=cfg.local.nTapers;
            end
            if strcmp(kind,'kalman')
                ap.pStayOn=[.98 .995 .999]; ap.pStayOff=[.98 .995 .999];
                ap.processFraction=[.0003 .001 .003]; ap.studentNu=[3 5 10];
                ap.observationFraction=[.05 .15 .30]; ap.damping=[.995 .999 .9995];
                defaults.pStayOn=.995; defaults.pStayOff=.995; defaults.processFraction=.001;
                defaults.studentNu=5; defaults.observationFraction=.15; defaults.damping=.999;
            end
        end
    case 'spectrum'
        pp.windowSeconds=[2 3 4 6]; pp.spectralOverlap=[0 .5 .75]; background=true;
        if strcmp(kind,'periodogram'), pp=rmfield(pp,{'windowSeconds','spectralOverlap'}); end
        if strcmp(kind,'long_welch'), pp.windowSeconds=[6 8 10]; end
        if any(strcmp(kind,{'multitaper','thomson'}))
            ap.NW=[1.5 2 2.5 3]; ap.nTapers=[1 2 3 4]; defaults.NW=2; defaults.nTapers=3;
        end
        if strcmp(kind,'thomson')
            threshold=false; background=false; ap.lineFThreshold=[4 8 12 20]; defaults.lineFThreshold=10;
        elseif strcmp(kind,'trimmed')
            ap.trimFraction=[.05 .10 .20]; defaults.trimFraction=.10;
        elseif strcmp(kind,'irasa')
            background=false;
            ap.resamplingFactors={ [1.1 1.3 1.5], [1.1 1.3 1.5 1.7 1.9] };
            defaults.resamplingFactors=[1.1 1.3 1.5 1.7 1.9];
        elseif any(strcmp(kind,{'acf','cepstrum','lomb'}))
            threshold=false; background=false; ap.alternativeThreshold=[2 4 6 10]; defaults.alternativeThreshold=4;
            if strcmp(kind,'acf'), ap.acfReference=[.05 .10 .20]; defaults.acfReference=.10; end
            if strcmp(kind,'lomb'), ap.frequencyStep=[.025 .05 .10]; defaults.frequencyStep=.05; end
        elseif any(strcmp(kind,{'robust_background','filter_aware','specparam'}))
            background=false;
        end
    case 'ar'
        ap.order=[4 8 12 20]; defaults.order=8;
        if strcmp(kind,'arma'), ap.maOrder=[1 3 5]; defaults.maOrder=3; end
        if strcmp(kind,'tvar')
            threshold=false; ap.forgetting=[.98 .995 .999]; ap.predictionRatio=[.10 .20 .40];
            defaults.forgetting=.995; defaults.predictionRatio=.20; ap.minCycles=[3 4 5];
        else, background=true;
        end
    case 'subspace'
        threshold=false; ap.embedding=[24 40 64 96]; ap.signalRank=[2 4 6];
        ap.lineFThreshold=[4 6 10 16]; defaults.embedding=40; defaults.signalRank=2; defaults.lineFThreshold=6;
        if strcmp(kind,'atomic')
            ap=rmfield(ap,{'embedding','signalRank'});
            ap.atomicSamples=[128 192 256]; ap.atomicPenalty=[.10 .30 .50 .80];
            ap.iterations=[40 80 160]; ap.rho=[.5 1 2]; ap.atomicAggregation={'mean','maximum'};
            defaults.atomicSamples=192; defaults.atomicPenalty=.5; defaults.iterations=80;
            defaults.rho=1; defaults.atomicAggregation='maximum';
        end
    case 'timefrequency'
        local=true; pp.localHop=[.1 .2 .3]; ap.minCycles=[3 4 5];
        ap.switchPenalty=[.5 1 2]; ap.frequencyPenalty=[.15 .35 .7];
        if any(strcmp(kind,{'multires','viterbi'}))
            ap.longWindow=[2 3 4]; defaults.longWindow=3; background=true;
            if strcmp(kind,'multires'), ap.shortWindow=[.75 1 1.5]; defaults.shortWindow=1; end
        else
            ap.frequencyStep=[.125 .25 .5]; defaults.frequencyStep=.25;
            if any(strcmp(kind,{'morlet','morse','bosc','ebosc','fbosc','obosc'}))
                ap.cycles=[3 4 5 7]; defaults.cycles=5;
            end
            if any(strcmp(kind,{'morse','bosc','ebosc','fbosc','obosc'}))
                ap.morseGamma=[2 3 4]; defaults.morseGamma=3;
            end
            if any(strcmp(kind,{'superlet','stppto'}))
                ap.cycleCounts={ [2 4 6], [3 6 9], [4 8 12] }; defaults.cycleCounts=[3 6 9];
            end
            if any(strcmp(kind,{'sst','reassign','conceft'}))
                ap.windowSigma=[.15 .25 .40]; defaults.windowSigma=.25;
            end
            if any(strcmp(kind,{'sst','conceft'}))
                ap.coefficientFloor=[.10 .15 .30]; defaults.coefficientFloor=.15;
            end
            if strcmp(kind,'conceft'), ap.projections=[4 8 16]; defaults.projections=8; end
            if strcmp(kind,'stransform'), ap.stockwellScale=[.75 1 1.5]; defaults.stockwellScale=1; end
            if strcmp(kind,'spwvd'), ap.lagSeconds=[.15 .30 .50]; defaults.lagSeconds=.30; end
            if strcmp(kind,'ebosc'), ap.edgeCycles=[0 .25 .5]; defaults.edgeCycles=.5; end
            background=~any(strcmp(kind,{'bosc','ebosc','fbosc','obosc','stppto'}));
        end
    case 'representation'
        local=true; pp.localWindowSeconds=[1 1.5 2 3]; pp.localHop=[.1 .2 .3];
        ap.minCycles=[3 4 5]; ap.switchPenalty=[.5 1 2]; ap.frequencyPenalty=[.15 .35 .7]; background=true;
        if strcmp(kind,'robust'), ap.winsorFraction=[.01 .02 .05]; defaults.winsorFraction=.02; end
    case 'cyclic'
        ap.segmentSeconds=[.25 .5 1]; ap.segmentHop=[.125 .25 .5];
        ap.coherenceFloor=[.0025 .005 .01]; ap.minimumCoherence=[.01 .025 .05];
        defaults.segmentSeconds=.5; defaults.segmentHop=.25; defaults.coherenceFloor=.005; defaults.minimumCoherence=.025;
        if ~strcmp(kind,'global')
            local=true; pp.localHop=[.1 .2 .3]; ap.cyclicWindow=[1.5 3 5]; defaults.cyclicWindow=3;
            ap.minCycles=[3 4 5]; ap.switchPenalty=[.5 1 2]; ap.frequencyPenalty=[.15 .35 .7];
        end
    case 'auxiliary'
        if strcmp(kind,'bicoherence')
            threshold=false; pp.windowSeconds=[2 3 4 6]; ap.bicoherenceThreshold=[.1 .2 .4];
            defaults.bicoherenceThreshold=.2; background=true;
        else
            local=true; pp.localWindowSeconds=[1 1.5 2 3]; pp.localHop=[.1 .2 .3];
            ap.minCycles=[3 4 5]; ap.switchPenalty=[.5 1 2]; ap.frequencyPenalty=[.15 .35 .7]; background=true;
        end
    otherwise
        error('tbench:UnsupportedSpace','No active parameter contract for %s.',family);
end
if threshold, ap.threshold=[2 3 4 6 8]; defaults.threshold=cfg.cfar.ratio; end
if local
    defaults.minCycles=cfg.local.minCycles; defaults.switchPenalty=cfg.local.switchPenalty;
    defaults.frequencyPenalty=cfg.local.frequencyPenalty;
end
if background
    ap.cfarGuardHz=[.5 .75 1]; ap.cfarReferenceHz=[3 5 7]; ap.cfarQuantile=[.4 .5 .7];
    defaults.cfarGuardHz=cfg.cfar.guardHz; defaults.cfarReferenceHz=cfg.cfar.referenceHz;
    defaults.cfarQuantile=cfg.cfar.quantile;
end
if tbench.option(p.preprocessParams,'freezeTransform',false)
    for name={'NW','nTapers','cycles','morseGamma','morseBeta','frequencyStep','cycleCounts', ...
            'windowSigma','coefficientFloor','projections','lagSeconds','stockwellScale','longWindow','shortWindow'}
        key=name{1}; if isfield(ap,key), ap=rmfield(ap,key); end
    end
    description.transformPolicy='frozen_preprocessor_from_explicit_cross_experiment';
else, description.transformPolicy='matched_preprocessor_transform_parameters_jointly_tuned';
end
for k=1:numel(cfg.tuning.spaces)
    s=cfg.tuning.spaces(k);
    if strcmp(s.analysisId,p.analysisId)
        if isfield(s,'preprocess')
            if ~isstruct(s.preprocess)||~all(ismember(fieldnames(s.preprocess),fieldnames(p.preprocessParams)))
                error('tbench:SearchSpace','Unknown preprocessing parameter in %s.',p.analysisId);
            end
            pp=s.preprocess;
        end
        if isfield(s,'analysis')
            allowed=union(fieldnames(ap),union(fieldnames(defaults),fieldnames(p.analysisParams)));
            if ~isstruct(s.analysis)||~all(ismember(fieldnames(s.analysis),allowed))
                error('tbench:SearchSpace','Unknown or unconnected analysis parameter in %s. Declare its implementation and space explicitly.',p.analysisId);
            end
            ap=s.analysis;
        end
    end
end
if tbench.option(p.preprocessParams,'freezeTransform',false)&& ...
        any(ismember(fieldnames(ap),{'NW','nTapers','cycles','morseGamma','morseBeta','frequencyStep','cycleCounts', ...
        'windowSigma','coefficientFloor','projections','lagSeconds','stockwellScale','longWindow','shortWindow'}))
    error('tbench:FrozenTransform','A cross experiment freezes its preprocessor transform. Change the preprocessor recipe instead of tuning ignored transform parameters.');
end
if isfield(ap,'ridge'), error('tbench:SearchSpace','Use tuning.ridgeValues for the clinical readout.'); end
end
