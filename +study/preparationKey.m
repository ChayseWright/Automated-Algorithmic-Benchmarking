function value = preparationKey(p)
%PREPARATIONKEY Analysis thresholds do not invalidate an identical transform.
params=p.analysisParams;
transformKeys={'NW','nTapers','cycles','morseBeta','morseGamma','frequencyStep','cycleCounts', ...
    'windowSigma','coefficientFloor','projections','lagSeconds','stockwellScale','longWindow','shortWindow'};
extra=setdiff(fieldnames(params),transformKeys); if ~isempty(extra), params=rmfield(params,extra); end
variant=p.preprocessAnalysisId;
if tbench.option(p.preprocessParams,'freezeTransform',false), params=p.preprocessParams.transformParams; end
if strcmp(p.inputKind,'local_spectrum')&&~strcmp(p.preprocessId,'multiresolution_spectra'), variant='shared_local_spectrum'; end
if any(strcmp(p.inputKind,{'carrier_time','envelope_time','positive_energy'})), variant='shared_time'; params=struct(); end
value=tbench.study.hash(struct('recipe',p.preprocessId,'implementation',p.implementationHash, ...
    'pp',p.preprocessParams,'transform',params,'variant',variant));
end
