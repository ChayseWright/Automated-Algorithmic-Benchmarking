function output = START_HERE(action,dataFolder,clinicalCSV,outputFolder,motionFolder)
%START_HERE Validate, demonstrate, inventory, or run the study.
% START_HERE('validate')
% START_HERE('demo')
% START_HERE('inventory','/data/MAT','/data/clinical.csv')
% study=START_HERE('run','/data/MAT','/data/clinical.csv','/results')
if nargin<1, action='validate'; end
root=fileparts(mfilename('fullpath')); addpath(root); addpath(fullfile(root,'tests'));
switch lower(char(action))
    case 'validate'
        output=run_native_validation(true);
    case 'demo'
        run(fullfile(root,'examples','run_study_demo.m'));
        output=study;
    case 'motion_demo'
        run(fullfile(root,'examples','run_motion_demo.m')); output=study;
    case {'run','inventory','motion'}
        if nargin<3, error('tbench:Usage','Supply the MAT-data folder and clinical CSV path.'); end
        cfg=tbench.studyConfig();
        if nargin>=4&&~isempty(outputFolder)
            cfg.outputDir=outputFolder; cfg.cacheDir=fullfile(outputFolder,'cache');
        end
        % Set timeScale, muscleNames, severityMode and synchronized here ONLY
        % according to your acquisition/clinical metadata; see README.md.
        if strcmpi(action,'inventory')
            subjects=tbench.attachClinical(tbench.subjectsFromFolder(dataFolder),clinicalCSV);
            output=tbench.inventory(subjects,cfg);
        elseif strcmpi(action,'motion')
            if nargin<5, motionFolder=''; end
            output=run_motion_benchmark(dataFolder,clinicalCSV,motionFolder,cfg);
        else, output=run_benchmark(dataFolder,clinicalCSV,cfg);
        end
    otherwise, error('tbench:Usage','Action must be validate, demo, motion_demo, inventory, run or motion.');
end
end
