function value = codeHash()
%CODEHASH Include executable source and the module manifest in cache identity.
root=fileparts(fileparts(fileparts(mfilename('fullpath'))));
files=[dir(fullfile(root,'+tbench','**','*.m'));dir(fullfile(root,'methods_manifest.json'))];
names=sort({files.name}); parts=cell(numel(files),1); %#ok<NASGU>
paths=arrayfun(@(f)fullfile(f.folder,f.name),files,'UniformOutput',false);
[paths,order]=sort(paths); files=files(order); %#ok<NASGU>
for k=1:numel(paths)
    relative=strrep(strrep(paths{k},root,''),'\','/');
    parts{k}=struct('path',relative,'sha',tbench.study.fileHash(paths{k}));
end
value=tbench.study.hash(parts);
end
