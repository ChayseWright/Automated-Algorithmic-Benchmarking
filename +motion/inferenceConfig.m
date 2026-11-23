function out = inferenceConfig(cfg)
%INFERENCECONFIG Do not expose motion files, references or clinical labels.
out=cfg;
for name={'motion','motionRuntime','motionProtocol','clinical','tuning'}
    key=name{1}; if isfield(out,key), out=rmfield(out,key); end
end
end
