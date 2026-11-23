function export(study)
%EXPORT Flat clinical tables, full parameters and stable run provenance.
folder=study.outputDir;
if isfield(study,'motion')
    motionNames=fieldnames(study.motion);
    for k=1:numel(motionNames)
        key=motionNames{k}; tbench.writeCSV(fullfile(folder,['motion_' key '.csv']),labels(study.motion.(key),study));
    end
end
names={'installation','inventory','reduction','outerPlan','clinicalPredictions','subjectFeatures', ...
    'levelAvailability','optimizationTrace','clinicalScorecard','pairedContrasts','clinicalRetention','tuningProvenance', ...
    'calibrationTraining','calibrationTest','physiologyScores','physiologyContrasts','physiologyRetention', ...
    'transferAudit','rhythmSummary','qualitySummary','phenotypeScorecard','commonEligibleScorecard', ...
    'operatingCurves','severityCaptureCurves','withinFoldAUC'};
for k=1:numel(names), name=names{k}; tbench.writeCSV(fullfile(folder,[name '.csv']),labels(study.(name),study)); end
names=fieldnames(study.repeatability);
for k=1:numel(names), name=names{k}; tbench.writeCSV(fullfile(folder,['B6_' name '.csv']),labels(study.repeatability.(name),study)); end
names=fieldnames(study.explanatory);
for k=1:numel(names), name=names{k}; tbench.writeCSV(fullfile(folder,['B5_' name '_summary.csv']),labels(study.explanatory.(name),study)); end
scores=struct([]);
for k=1:numel(study.selectionScores)
    s=study.selectionScores(k); row=s.metrics; row.pipelineId=s.pipelineId; row.candidateId=s.candidateId;
    row.candidateHash=s.candidateHash; row.configurationHash=s.configurationHash; row.outerFold=s.outerFold; row.eligible=s.eligible; row.violation=s.violation;
    row.searchRound=s.searchRound; row.selectionScope=s.selectionScope;
    row.threshold=s.threshold; row.ridge=s.ridge; row.strongDetection=s.calibration.strongDetection;
    row.strongPowerError=s.calibration.strongPowerError; row.calibrationCoverage=s.calibration.coverage;
    row.strongFrequencyError=s.calibration.strongFrequencyError; row.nullAlarmFraction=s.calibration.nullAlarmFraction;
    row.strongCalibrationCases=s.calibration.strongCases; row.nullCalibrationCases=s.calibration.nullCases;
    row.weakCalibrationCases=s.calibration.weakCases; row.weakFrequencyRecovered=s.calibration.weakFrequencyRecovered;
    row.lateStrongCases=s.calibration.lateStrongCases; row.lateStrongRecovered=s.calibration.lateStrongRecovered;
    row.motionUtility=s.motion.gain; row.motionReadoutRidge=s.motion.ridge;
    row.motionValidationSubjects=s.motion.subjects; row.motionInputCoverage=s.motion.inputCoverage;
    if isempty(scores), scores=row; else, scores(end+1)=row; end %#ok<AGROW>
end
tbench.writeCSV(fullfile(folder,'tuning_scores.csv'),labels(scores,study));
catalog=struct([]);
for k=1:numel(study.candidateCatalog)
    p=study.candidateCatalog(k);
    row=struct('pipelineId',p.id,'analysisId',p.analysisId,'preprocessId',p.preprocessId,'inputKind',p.inputKind, ...
        'candidateId',p.candidateId,'candidateHash',p.hash, ...
        'preprocessParameters',jsonencode(p.preprocessParams),'analysisParameters',jsonencode(p.analysisParams));
    if isempty(catalog), catalog=row; else, catalog(end+1)=row; end %#ok<AGROW>
end
tbench.writeCSV(fullfile(folder,'candidate_catalog.csv'),labels(catalog,study));
[spaces,stability]=tbench.study.parameterTables(study);
tbench.writeCSV(fullfile(folder,'parameter_spaces.csv'),labels(spaces,study));
tbench.writeCSV(fullfile(folder,'parameter_stability.csv'),labels(stability,study));
C=tbench.methodCatalog(); tbench.writeCSV(fullfile(folder,'method_availability.csv'),labels(C,study));
selected=struct([]);
for k=1:numel(study.selection)
    w=study.selection(k); p=study.candidateCatalog(w.candidateIndex);
    r=struct('scope','outer_training_selection','pipelineId',p.id,'outerFold',w.outerFold, ...
        'candidateId',p.candidateId,'candidateHash',p.hash,'configurationHash',w.configurationHash, ...
        'ridge',w.ridge,'threshold',w.threshold,'selectionStatus',w.selectionStatus, ...
        'motionReadoutRidge',w.motion.ridge,'motionValidationUtility',w.motion.gain, ...
        'preprocessParameters',jsonencode(p.preprocessParams),'analysisParameters',jsonencode(p.analysisParams));
    if isempty(selected), selected=r; else, selected(end+1)=r; end %#ok<AGROW>
end
for k=1:numel(study.deployment)
    w=study.deployment(k); p=w.pipeline;
    r.scope='all_training_deployment'; r.pipelineId=p.id; r.outerFold=0;
    r.candidateId=p.candidateId; r.candidateHash=p.hash; r.configurationHash=w.configurationHash;
    r.ridge=w.ridge; r.threshold=w.threshold; r.selectionStatus=w.selectionStatus;
    r.motionReadoutRidge=w.motionPolicy.ridge; r.motionValidationUtility=w.motionPolicy.gain;
    r.preprocessParameters=jsonencode(p.preprocessParams); r.analysisParameters=jsonencode(p.analysisParams);
    selected(end+1)=r; %#ok<AGROW>
end
tbench.writeCSV(fullfile(folder,'selected_parameters.csv'),labels(selected,study));
for name={'features','events'}
    path=fullfile(folder,[name{1} '.csv']); if ~exist(path,'file'), tbench.writeCSV(path,struct([])); end
end
tbench.study.atomicSave(fullfile(folder,'study.mat'),study);
manifest=struct('runId',study.runId,'runNumber',study.runNumber,'status',study.status, ...
    'protocolVersion',study.protocolVersion,'codeHash',study.codeHash,'config',study.config, ...
    'manifestAudit',study.manifestAudit,'sourceHashes',{study.sourceHashes}, ...
    'selection',study.selection,'cacheHits',study.cacheHits,'cacheMisses',study.cacheMisses);
manifest.motionEnabled=study.config.motion.enabled;
manifest.motionReferenceHashes=study.config.motionRuntime.subjectHashes;
tbench.study.writeJSON(fullfile(folder,'run_manifest.json'),manifest);
tbench.study.report(study);
end
function rows=labels(rows,study)
for k=1:numel(rows), rows(k).runId=study.runId; rows(k).runNumber=study.runNumber; end
end
