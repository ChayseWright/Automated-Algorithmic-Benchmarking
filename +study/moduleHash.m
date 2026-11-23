function value = moduleHash(p)
%MODULEHASH Include named module source and declared external dependencies.
parts=cell(1,3+numel(p.dependencyFiles)); handles={p.preprocess,p.analyze,p.benchmark};
for k=1:numel(handles)
    if ~isa(handles{k},'function_handle'), error('tbench:ModuleHandle','Module hooks must be function handles.'); end
    info=functions(handles{k});
    if strcmp(info.type,'anonymous')||strcmp(info.type,'nested')
        error('tbench:OpaqueModule','Use a named function file and explicit pipeline parameters instead of captured module state.');
    end
    file=info.file;
    if isempty(file)||exist(file,'file')~=2, error('tbench:ModuleSource','Module source file is required: %s.',func2str(handles{k})); end
    parts{k}=struct('entrypoint',func2str(handles{k}),'sha256',tbench.study.fileHash(file));
end
for j=1:numel(p.dependencyFiles)
    [~,name,extension]=fileparts(p.dependencyFiles{j});
    parts{3+j}=struct('entrypoint',[name extension],'sha256',tbench.study.fileHash(p.dependencyFiles{j}));
end
value=tbench.study.hash(parts);
end
