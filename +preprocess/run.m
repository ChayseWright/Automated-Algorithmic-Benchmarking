function prepared = run(c,cfg,pipeline)
%RUN Matched input representation; reference y/envelope remain unchanged.
p=pipeline.preprocessParams;
owner=tbench.option(pipeline,'preprocessAnalysisId',pipeline.analysisId);
if isempty(owner), owner=pipeline.analysisId; end
transform=pipeline.analysisParams;
if tbench.option(p,'freezeTransform',false), transform=p.transformParams; end
transformKeys={'NW','nTapers','cycles','morseBeta','morseGamma','frequencyStep','cycleCounts', ...
    'windowSigma','coefficientFloor','projections','lagSeconds','stockwellScale','longWindow','shortWindow'};
extra=setdiff(fieldnames(transform),transformKeys);
if ~isempty(extra), transform=rmfield(transform,extra); end
cfg.localHop=p.localHop; cfg.spectralOverlap=p.spectralOverlap;
parameters=[p.rawHighpassHz p.envelopePassHz p.envelopeStopHz];
if isfield(c,'cachedDetector')&&isequal(c.cachedDetectorParameters,parameters)
    prepared=c.cachedDetector;
else, prepared=tbench.preprocess.detector(c,cfg,p);
end
prepared.inputKind=pipeline.inputKind; prepared.preprocessingId=pipeline.preprocessId;
detector=prepared;
if strcmp(pipeline.inputKind,'carrier_time')
    if isfield(c,'reducedRaw')
        prepared.raw=c.reducedRaw; prepared.rawTime=c.reducedRawTime; prepared.rawFs=c.reducedRawFs;
        prepared.rawValid=isfinite(prepared.raw);
    end
elseif strcmp(pipeline.inputKind,'envelope_spectrum')
    seconds=p.windowSeconds; taper=p.taper; input=detector.envelope;
    if strcmp(owner,'periodogram')
        R=tbench.validRuns(detector.valid&isfinite(input)); seconds=0;
        if ~isempty(R)
            [~,j]=max(R(:,2)-R(:,1)); input=input(R(j,1):R(j,2)); seconds=numel(input)/c.fs;
        end
    elseif strcmp(owner,'long_welch'), seconds=min(p.windowSeconds,numel(c.t)/c.fs);
    end
    if seconds>0
        NW=tbench.option(transform,'NW',2);
        K=tbench.option(transform,'nTapers',3);
        [P,f,Z]=tbench.spectrum(input,c.fs,seconds,cfg.spectralOverlap,taper,NW,K);
        prepared.preparedSpectrum=struct('P',P,'f',f,'Z',Z,'seconds',seconds,'taper',taper,'NW',NW,'nTapers',K);
    end
elseif strcmp(pipeline.inputKind,'local_spectrum')&&strcmp(pipeline.preprocessId,'multiresolution_spectra')
    [P,f,t,good]=tbench.localEvidence(detector,cfg,tbench.option(transform,'longWindow',3),'hann');
    prepared.preparedMultiEvidence=struct('long',struct('P',P,'f',f,'t',t,'good',good),'short',[]);
    if strcmp(owner,'multiresolution_stft')
        [P,f,t,good]=tbench.localEvidence(detector,cfg,tbench.option(transform,'shortWindow',1),'hann');
        prepared.preparedMultiEvidence.short=struct('P',P,'f',f,'t',t,'good',good);
    end
elseif strcmp(pipeline.inputKind,'local_spectrum')
    local=cfg; local.local.NW=tbench.option(transform,'NW',cfg.local.NW);
    local.local.nTapers=tbench.option(transform,'nTapers',cfg.local.nTapers);
    [P,f,t,good]=tbench.localEvidence(detector,local,p.localWindowSeconds,p.taper);
    prepared.preparedEvidence=struct('P',P,'f',f,'t',t,'good',good,'taper',p.taper);
elseif strcmp(pipeline.inputKind,'complex_tf')
    C=tbench.methodCatalog(); j=find(strcmp({C.id},owner),1); kind=C(j).variant;
    if any(strcmp(kind,{'bosc','ebosc','fbosc','obosc'})), kind='morse'; end
    [P,f,t,good,W]=tbench.tfMap(detector,cfg,kind,transform);
    prepared.preparedTF=struct('P',P,'f',f,'t',t,'good',good,'W',W,'transform',kind);
end
prepared.preprocessAudit=struct('id',pipeline.preprocessId,'inputKind',pipeline.inputKind, ...
    'referenceFs',c.fs,'detectorRawFs',prepared.rawFs,'parameters',p,'transformOwner',owner,'transformParameters',transform);
end
