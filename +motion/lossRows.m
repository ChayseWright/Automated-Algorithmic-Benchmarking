function [rows,models] = lossRows(view,trainingIds,testIds,cfg,level,lambda)
%LOSSROWS Subject-held-out mechanical-amplitude readout; no test-target fitting.
rows=struct([]); models=struct([]); if isempty(view)||isempty(view.rows), return; end
if ~isempty(intersect(trainingIds,testIds)), error('tbench:MotionLeakage','Motion readout training/test subjects overlap.'); end
A=view.rows; names=unique({A.target});
for k=1:numel(names)
    use=strcmp({A.target},names{k});
    if ~isempty(cfg.motion.primaryTargets)&&~any(strcmp(cfg.motion.primaryTargets,names{k})), continue; end
    target=A(use);
    if numel(unique({target.unit}))~=1||numel(unique({target.kind}))~=1||numel(unique({target.frame}))~=1
        error('tbench:MotionTargetConsistency','A named motion target changes coordinate kind or physical units.');
    end
    train=ismember({target.subject},trainingIds); test=ismember({target.subject},testIds);
    power=[target.meanPower]'; power(isfinite(power)&power<0)=0;
    y=sqrt(power); % Missing reference power must never become a zero target.
    eligible=isfinite(y)&[target.referenceValidFraction]'>=cfg.motion.minimumReferenceCoverage;
    if sum(train(:)&eligible)<cfg.motion.minimumFitRows||~any(test), continue; end
    if numel(unique({target(train&eligible').subject}))<cfg.motion.minimumMotionSubjects, continue; end
    taskX=taskDesign(target,cfg.tasks); base=vertcat(target.B0);
    rawLevel=level; frozen=strncmp(level,'fixed_',6);
    if frozen, rawLevel='R'; end
    X=vertcat(target.(rawLevel));
    inputOK=[target.inputCoverage]'>=cfg.motion.minimumSourceCoverage&all(isfinite(X),2);
    tr=train(:)&eligible; te=test(:)&eligible;
    yScale=positiveScale(y(tr)); z=asinh(y/yScale);
    xScale=scales(X(tr&inputOK,:)); bScale=scales(base(tr,:));
    transformed=asinh(bsxfun(@rdivide,X,xScale)); B=asinh(bsxfun(@rdivide,base,bScale));
    weights=subjectWeights(target,cfg.tasks);
    baseline=tbench.motion.fitLinear([taskX B],z,cfg.motion.baselineRidge,weights.*double(tr));
    model=tbench.motion.fitLinear([taskX transformed],z,lambda,weights.*double(tr&inputOK));
    if isempty(baseline), continue; end
    replacementValues=nan(1,size(X,2));
    if frozen
        short=level(7:end); retained=retainedColumns(short);
        for d=setdiff(1:size(X,2),retained)
            v=X(tr&inputOK,d); v=v(isfinite(v)); replacement=0;
            if ~isempty(v), replacement=median(v); end
            replacementValues(d)=replacement;
            transformed(:,d)=asinh(replacement/xScale(d));
        end
    end
    if nargout>1
        fitted=struct('target',target(1).target,'kind',target(1).kind,'unit',target(1).unit, ...
            'frame',target(1).frame,'level',level,'ridge',lambda,'baselineRidge',cfg.motion.baselineRidge,'yScale',yScale, ...
            'xScale',xScale,'baselineScale',bScale,'model',model,'baseline',baseline, ...
            'replacementValues',replacementValues,'taskNames',{cfg.tasks}, ...
            'trainingSubjects',{trainingIds},'targetQuantity','guarded_band_RMS_amplitude');
        models=append(models,fitted);
    end
    predictionB=tbench.motion.predictLinear(baseline,[taskX B]);
    null=mean(z(tr)); predictionB(~isfinite(predictionB))=null;
    prediction=tbench.motion.predictLinear(model,[taskX transformed]);
    prediction(~inputOK|~isfinite(prediction))=predictionB(~inputOK|~isfinite(prediction));
    % Nonnegative physical amplitudes; the same bound applies to both readouts.
    prediction=max(0,prediction); predictionB=max(0,predictionB);
    physical=yScale*sinh(min(50,prediction)); physicalB=yScale*sinh(min(50,predictionB));
    variance=max(cfg.motion.modelVarianceFloor,baseline.variance);
    for j=find(te(:))'
        r=struct('subject',target(j).subject,'task',target(j).task,'trial',target(j).trial, ...
            'target',target(j).target,'kind',target(j).kind,'unit',target(j).unit, ...
            'level',level,'ridge',lambda,'scope','unseen_subject_trial_RMS_amplitude', ...
            'referenceValue',y(j),'predictedValue',physical(j),'baselineValue',physicalB(j), ...
            'referenceTransformed',z(j),'predictionTransformed',prediction(j),'trainingAmplitudeScale',yScale, ...
            'physicalSquaredError',(y(j)-physical(j))^2,'baselinePhysicalSquaredError',(y(j)-physicalB(j))^2, ...
            'transformedLoss',(z(j)-prediction(j))^2,'baselineTransformedLoss',(z(j)-predictionB(j))^2, ...
            'logScoreGain',((z(j)-predictionB(j))^2-(z(j)-prediction(j))^2)/(2*variance), ...
            'inputAvailable',inputOK(j),'referenceCoverage',target(j).referenceValidFraction, ...
            'trainingSubjects',strjoin(trainingIds,' '),'referenceHash',target(j).referenceHash);
        rows=append(rows,r);
    end
end
end
function A=taskDesign(rows,tasks)
A=zeros(numel(rows),max(0,numel(tasks)-1));
for j=2:numel(tasks), A(:,j-1)=strcmp({rows.task},tasks{j})'; end
end
function s=positiveScale(x)
x=x(isfinite(x)&x>0); s=1; if ~isempty(x), s=median(x); end
end
function s=scales(X)
s=ones(1,size(X,2)); for j=1:size(X,2), s(j)=positiveScale(abs(X(:,j))); end
end
function w=subjectWeights(rows,tasks)
w=zeros(numel(rows),1); ids=unique({rows.subject});
for s=1:numel(ids)
    for type='PK'
        names=tasks(cellfun(@(x)x(1)==type,tasks));
        for j=1:numel(names)
            ix=strcmp({rows.subject},ids{s})&strcmp({rows.task},names{j});
            if any(ix), w(ix)=.5/(numel(names)*sum(ix)); end
        end
    end
end
end
function retained=retainedColumns(level)
switch level
    case 'F0', retained=[15 3];
    case 'F1', retained=1:3;
    case 'F2', retained=1:6;
    case 'F3', retained=1:9;
    otherwise, error('tbench:MotionCompression','Unknown frozen feature level.');
end
end
function a=append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
