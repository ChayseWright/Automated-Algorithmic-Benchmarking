function result = run_native_validation(includeEndToEnd)
%RUN_NATIVE_VALIDATION Native parser audit plus legacy, study and motion tests.
if nargin<1, includeEndToEnd=true; end
root=fileparts(fileparts(mfilename('fullpath'))); addpath(root); addpath(fullfile(root,'tests'));
files=dir(fullfile(root,'**','*.m')); issues=struct('file',{},'line',{},'message',{});
for k=1:numel(files)
    file=fullfile(files(k).folder,files(k).name); M=checkcode(file,'-id');
    for j=1:numel(M)
        issues(end+1)=struct('file',file,'line',M(j).line,'message',M(j).message); %#ok<AGROW>
    end
end
if ~isempty(issues)
    bad=~cellfun('isempty',regexpi({issues.message},'parse error|syntax error|invalid expression|illegal use of reserved keyword'));
    assert(~any(bad),'tbench:NativeSyntax','Native checkcode reported a syntax/parse error.');
end
result=struct('engine',version,'time',datestr(now,30),'checkcodeIssues',issues, ...
    'legacy',run_selftests(),'study',run_study_tests(includeEndToEnd,false), ...
    'motion',run_motion_tests(includeEndToEnd));
fprintf('Native validation completed in %s. No patient-performance claim follows from synthetic tests.\n',version);
end
