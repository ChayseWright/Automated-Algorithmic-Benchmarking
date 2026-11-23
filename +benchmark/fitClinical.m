function model = fitClinical(X,subjects,available,cfg,lambda)
%FITCLINICAL Clinical labels enter only this cohort readout.
y=strcmp({subjects.group},'ET')'; valid=available(:)&all(isfinite(X),2);
model.featureTransform=powerTransform(X,available(:),~y,cfg.clinical.featureLevel);
X=applyPowerTransform(X,model.featureTransform);
if sum(valid)<2||numel(unique(y(valid)))<2
    % Still provide a baseline prediction, but no eligible fitted classifier.
    valid=true(size(y)); X(:)=0; fitStatus='baseline_only';
else, fitStatus='ok';
end
model.diagnosis=tbench.benchmark.fitLogistic(X(valid,:),y(valid),lambda,true);
model.status=fitStatus; model.trainingSubjects={subjects(valid).id};
model.trainingPrior=.5; model.lambda=lambda; % Balanced working prior, not clinic prevalence.
[group,cuts,label]=tbench.benchmark.severityGroups(subjects,cfg);
model.severityCuts=cuts; model.severityDescription=label;
use=available(:)&y&group(:)>0&all(isfinite(X),2);
model.cumulative=cell(1,2);
for k=1:2
    if sum(use)>=2
        model.cumulative{k}=tbench.benchmark.fitLogistic(X(use,:),double(group(use)<=k),lambda,false);
    end
end
% A numeric score readout is secondary; ordinal probabilities remain separate.
score=[subjects.severity_score]'; use=available(:)&y&isfinite(score)&all(isfinite(X),2);
model.numericSeverity=[];
if sum(use)>=3
    A=X(use,:); mu=mean(A,1); sd=std(A,1,1); sd(sd<1e-10)=1;
    B=[ones(sum(use),1) bsxfun(@rdivide,bsxfun(@minus,A,mu),sd)];
    penalty=lambda*eye(size(B,2)); penalty(1,1)=0;
    model.numericSeverity=struct('mu',mu,'sd',sd,'beta',(B'*B+penalty)\(B'*score(use)), ...
        'minimum',min(score(use)),'maximum',max(score(use)));
end
end
function transform=powerTransform(X,available,control,level)
switch level
    case 'B0', columns=[1 2];
    case 'F0', columns=1;
    case 'F1', columns=[1 2];
    case 'F2', columns=[1 2 5];
    case 'F3', columns=[1 2 5 7];
    otherwise, columns=[1 2 5 7 10 11 15];
end
columns=columns(columns<=size(X,2)); scale=ones(size(columns));
for k=1:numel(columns)
    v=X(available&control,columns(k)); v=v(isfinite(v)&v>0);
    if isempty(v), v=X(available,columns(k)); v=v(isfinite(v)&v>0); end
    if ~isempty(v), scale(k)=median(v); end
end
transform=struct('columns',columns,'scales',scale,'rule','training_control_positive_median_fallback_training');
end
function X=applyPowerTransform(X,transform)
for k=1:numel(transform.columns)
    j=transform.columns(k); X(:,j)=asinh(X(:,j)/transform.scales(k));
end
end
