function [subjects,audit] = preflight(subjects,cfg)
%PREFLIGHT Fail before expensive work when labels, modules or contracts conflict.
if isempty(subjects)||~isstruct(subjects), error('tbench:Subjects','Provide a nonempty subject manifest.'); end
if ~all(isfield(subjects,'id')), error('tbench:SubjectID','Every participant needs an explicit ID.'); end
if ~cfg.clinical.enabled, error('tbench:ClinicalRequired','Use the legacy runner for a study without clinical/control labels.'); end
tbench.motion.validateConfig(cfg.motion,cfg);
folds=[cfg.outerFolds cfg.tuning.innerFolds];
if any(~isfinite(folds)|folds<2|folds~=floor(folds)), error('tbench:Folds','Use integer inner/outer fold counts of at least two.'); end
if ~any(strcmp(cfg.clinical.featureLevel,{'B0','F0','F1','F2','F3','R'})), error('tbench:FeatureLevel','Unknown primary feature level.'); end
if numel(cfg.muscleNames)~=15, error('tbench:Muscles','Supply exactly 15 acquisition-order muscle names.'); end
if numel(unique(cfg.tasks))~=numel(cfg.tasks)||numel(unique(cfg.trials))~=numel(cfg.trials)|| ...
        ~any(cellfun(@(x)x(1)=='P',cfg.tasks))||~any(cellfun(@(x)x(1)=='K',cfg.tasks))
    error('tbench:TaskBalance','Use unique tasks/trials and both postural and kinetic tasks for this balanced study readout.');
end
if ~isequal(size(cfg.clinical.costMatrix),[4 4])||any(~isfinite(cfg.clinical.costMatrix(:)))||any(cfg.clinical.costMatrix(:)<0)
    error('tbench:ClinicalCosts','Supply a finite nonnegative 4-by-4 decision cost matrix.');
end
if numel(cfg.rawBand)~=2||any(~isfinite(cfg.rawBand))||cfg.rawBand(1)<0||cfg.rawBand(2)<=cfg.rawBand(1)
    error('tbench:CarrierBand','Provide two finite ordered carrier-band edges.');
end
fractions=[cfg.clinical.minSubjectCoverage cfg.clinical.minClinicalCoverage cfg.clinical.minSevereCapture ...
    cfg.tuning.minStrongDetection cfg.tuning.minCalibrationCoverage cfg.tuning.maxNullAlarmFraction];
if any(~isfinite(fractions)|fractions<0|fractions>1)||~isfinite(cfg.clinical.controlFPR)
    error('tbench:Fractions','Coverage and detection constraints must be finite fractions in [0,1].');
end
if numel(unique({cfg.pipelines.id}))~=numel(cfg.pipelines), error('tbench:PipelineIDs','Pipeline IDs must be unique.'); end
if any(cellfun('isempty',regexp({cfg.pipelines.id},'^[A-Za-z0-9_][A-Za-z0-9_-]*$','once')))
    error('tbench:PipelineIDs','Pipeline IDs may contain letters, digits, underscores and hyphens.');
end
if any(cellfun('isempty',regexp(cfg.tasks,'^[PK][1-7]$','once')))|| ...
        any(cellfun('isempty',regexp(cfg.trials,'^trial[1-3]$','once')))
    error('tbench:TaskNames','Use acquisition task names P1--P7/K1--K7 and trial1--trial3.');
end
if ~any(strcmp(cfg.clinical.priority,{'mild','severe'}))||~any(strcmp(cfg.clinical.powerDomain,{'mvc','absolute'}))
    error('tbench:ClinicalConfig','Check the clinical priority and power-domain names.');
end
if ~isfinite(cfg.clinical.upperRecordQuantile)||cfg.clinical.upperRecordQuantile<=0||cfg.clinical.upperRecordQuantile>1
    error('tbench:UpperBurden','Prespecify upperRecordQuantile in (0,1].');
end
if isempty(cfg.tuning.ridgeValues)||any(cfg.tuning.ridgeValues<=0|~isfinite(cfg.tuning.ridgeValues))|| ...
        cfg.tuning.maxCandidates<1||cfg.tuning.maxCandidates~=floor(cfg.tuning.maxCandidates)
    error('tbench:TuningConfig','Provide positive finite ridge values and a positive integer candidate budget.');
end
counts=[cfg.tuning.refinementRounds cfg.tuning.refinementCandidates cfg.tuning.refinementSeeds cfg.tuning.maxGridCandidates];
if any(~isfinite(counts)|counts~=floor(counts)|counts<0)||any(counts(2:4)<1)|| ...
        ~any(strcmp(cfg.tuning.searchMode,{'coverage','grid'}))
    error('tbench:SearchConfig','Use nonnegative integer refinement rounds, positive integer budgets/seeds, and coverage or grid search.');
end
if ~isscalar(cfg.execution.candidateBatchSize)||~isfinite(cfg.execution.candidateBatchSize)|| ...
        cfg.execution.candidateBatchSize<1||cfg.execution.candidateBatchSize~=floor(cfg.execution.candidateBatchSize)
    error('tbench:BatchSize','Use a positive integer candidateBatchSize.');
end
if cfg.bench.calibration
    patterns={'continuous','early','middle','late','multiple'};
    if cfg.calibration.nCases<1||cfg.calibration.nCases~=floor(cfg.calibration.nCases)|| ...
            isempty(cfg.calibration.patterns)||~all(ismember(cfg.calibration.patterns,patterns))|| ...
            any(cfg.calibration.depths<0|cfg.calibration.depths>=1)
        error('tbench:CalibrationConfig','Use a positive calibration count, supported burst patterns and modulation depths in [0,1).');
    end
end
ids=cellfun(@char,{subjects.id},'UniformOutput',false);
safe=cellfun(@(x)lower(regexprep(x,'[^A-Za-z0-9_-]','_')),ids,'UniformOutput',false);
if any(cellfun('isempty',ids))||numel(unique(safe))~=numel(ids), error('tbench:IDs','Subject IDs must have unique safe filenames.'); end
groups=cell(numel(subjects),1);
for k=1:numel(subjects)
    subjects(k).id=ids{k};
    group=lower(char(tbench.option(subjects(k),'group','')));
    if any(strcmp(group,{'control','healthy'})), group='control';
    elseif strcmp(group,'et'), group='ET';
    else, error('tbench:Group','Explicit ET/control group is required for %s.',ids{k});
    end
    subjects(k).group=group; groups{k}=group;
    value=tbench.option(subjects(k),'severity_score',NaN);
    if ~isnumeric(value)||~isreal(value)||~isscalar(value)||isinf(value), error('tbench:Severity','Real numeric severity_score or NaN is required.'); end
    subjects(k).severity_score=value;
    subjects(k).severity_kind=char(tbench.option(subjects(k),'severity_kind',cfg.clinical.severityKind));
    if isempty(subjects(k).severity_kind), subjects(k).severity_kind=cfg.clinical.severityKind; end
    subjects(k).severity_group=lower(char(tbench.option(subjects(k),'severity_group','')));
    subjects(k).phenotype=char(tbench.option(subjects(k),'phenotype','unknown'));
    if strcmp(group,'ET')&&strcmp(cfg.clinical.severityMode,'provided')&& ...
            ~any(strcmp(subjects(k).severity_group,{'mild','moderate','severe'}))
        error('tbench:SeverityClass','Provided severity mode needs mild/moderate/severe for every ET subject.');
    end
end
nc=sum(strcmp(groups,'control')); ne=sum(strcmp(groups,'ET'));
if numel(subjects)~=cfg.expectedSubjects
    warning('tbench:SubjectCount','Expected %d people; running the supplied %d and retaining actual denominators.',cfg.expectedSubjects,numel(subjects));
end
if cfg.clinical.enabled&&(nc<cfg.outerFolds||ne<cfg.outerFolds)
    error('tbench:SplitSize','At least outerFolds controls and ET subjects are required.');
end
if min([nc-ceil(nc/cfg.outerFolds) ne-ceil(ne/cfg.outerFolds)])<2
    error('tbench:InnerSize','Every outer-training set needs at least two controls and two ET participants.');
end
if ~strcmp(cfg.clinical.severityMode,'provided')&&~any(isfinite([subjects(strcmp(groups,'ET')).severity_score]))
    error('tbench:MissingSeverity','Join the numeric TETRAS target, or supply independently defined clinical severity groups.');
end
if ~cfg.bench.calibration&&cfg.tuning.calibrationGate
    error('tbench:CalibrationGate','Calibration is disabled but its eligibility gate is enabled.');
end
if cfg.bench.calibration&&(cfg.calibration.fs<=2*cfg.rawBand(2)||cfg.calibration.fs<cfg.analysisFs)
    error('tbench:CalibrationFs','Synthetic Fs must exceed twice the carrier high cutoff and be no lower than analysisFs.');
end
if any(strcmp(cfg.clinical.featureLevel,{'F2','F3'}))&&any(~[cfg.pipelines.expectedDynamic])
    error('tbench:PrimaryCapacity','An episodic primary feature level requires dynamic output from every compared pipeline.');
end
if cfg.bench.controlInjection&&any(cfg.injection.depths<=0|cfg.injection.depths>=1)
    error('tbench:InjectionDepth','Counterfactual modulation depth must lie strictly between 0 and 1.');
end
if cfg.analysisFs/2<=cfg.reduction.envelopeStopHz|| ...
        cfg.reduction.envelopePassHz<cfg.harmonics*cfg.searchBand(2)
    error('tbench:HarmonicBandwidth','Reporting passband must include all declared harmonics with an output Nyquist margin.');
end
if cfg.searchBand(1)<=cfg.baselineCutoff||cfg.clinical.controlFPR<0||cfg.clinical.controlFPR>=1
    error('tbench:Configuration','Check search-band and control FPR settings.');
end
if isempty(cfg.pipelines), error('tbench:Pipelines','Register at least one complete pipeline.'); end
for k=1:numel(cfg.pipelines)
    p=cfg.pipelines(k);
    if ~strcmp(p.readiness,'implemented')
        error('tbench:UnavailableModule','%s is %s; it cannot enter an automated study.',p.analysisId,p.readiness);
    end
    if isempty(which(['tbench.methods.' p.analysisId]))
        error('tbench:MissingModule','Missing public analysis entry %s.',p.analysisId);
    end
    tbench.study.validatePair(p);
    tbench.study.moduleHash(p);
end
audit=struct('subjects',numel(subjects),'ET',ne,'controls',nc, ...
    'plannedTrials',numel(subjects)*numel(cfg.tasks)*numel(cfg.trials), ...
    'plannedMuscleRecords',numel(subjects)*numel(cfg.tasks)*numel(cfg.trials)*15, ...
    'severityMode',cfg.clinical.severityMode,'priority',cfg.clinical.priority);
end
