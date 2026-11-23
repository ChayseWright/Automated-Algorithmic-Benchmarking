function [folder,runId,runNumber] = newRun(root,label)
%NEWRUN Number plus unique suffix prevents simultaneous-run collisions.
if ~exist(root,'dir'), mkdir(root); end
lockfile=fullfile(root,'.run_number_lock');
lock=javaObject('java.io.File',lockfile); begin=tic;
while ~lock.createNewFile()
    if toc(begin)>10, error('tbench:RunLock','Run-number lock is busy; check for a stale lock after another process ended.'); end
    pause(.05);
end
cleanup=onCleanup(@()delete(lockfile)); %#ok<NASGU>
files=dir(fullfile(root,'run_*')); numbers=0;
for k=1:numel(files)
    v=regexp(files(k).name,'^run_(\d+)_','tokens','once');
    if ~isempty(v), numbers(end+1)=str2double(v{1}); end %#ok<AGROW>
end
runNumber=max(numbers)+1; [~,nonce]=fileparts(tempname(root));
label=regexprep(char(label),'[^A-Za-z0-9_-]','_');
runId=sprintf('run_%04d_%s_%s_%s',runNumber,datestr(now,'yyyymmddTHHMMSS'),label,nonce);
folder=fullfile(root,runId); mkdir(folder);
end
