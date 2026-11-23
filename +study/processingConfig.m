function out = processingConfig(cfg)
%PROCESSINGCONFIG Strip clinical labels/readouts and output locations from keys.
out=cfg;
remove={'outputDir','cacheDir','clinical','tuning','pipelines','execution', ...
    'methods','parameterSelector','bootstrapSamples','expectedSubjects','verbose','workDir','runId','runNumber', ...
    'bench','calibration','injection','temporal','prediction','motion','motionRuntime','motionProtocol'};
for k=1:numel(remove), if isfield(out,remove{k}), out=rmfield(out,remove{k}); end; end
end
