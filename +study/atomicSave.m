function atomicSave(filename,payload)
%ATOMICSAVE Write a complete MAT before publishing its cache/checkpoint name.
folder=fileparts(filename); if ~exist(folder,'dir'), mkdir(folder); end
temporary=[tempname(folder) '.mat']; cleanup=onCleanup(@()removeTemp(temporary)); %#ok<NASGU>
save(temporary,'payload','-v7');
[ok,message]=movefile(temporary,filename,'f'); if ~ok, error('tbench:Save','%s',message); end
end
function removeTemp(path)
if exist(path,'file'), delete(path); end
end
