function rows = calibrationCached(candidates,cfg,codeHash,offset)
%CALIBRATIONCACHED Shared by initial screening and training-only refinement.
if ~cfg.bench.calibration||isempty(candidates), rows=struct([]); return; end
key=tbench.study.hash(struct('code',codeHash,'candidates',{arrayfun(@(p)p.hash,candidates,'UniformOutput',false)}, ...
    'processing',tbench.study.processingConfig(cfg),'calibration',cfg.calibration,'offset',offset,'seed',cfg.seed));
filename=fullfile(cfg.cacheDir,'calibration',[key '.mat']);
if cfg.execution.cache&&exist(filename,'file')
    try
        S=load(filename,'payload');
        if strcmp(S.payload.key,key), rows=S.payload.rows; return; end
    catch
    end
end
rows=tbench.study.calibrate(candidates,cfg,offset);
if cfg.execution.cache, tbench.study.atomicSave(filename,struct('key',key,'rows',rows)); end
end
