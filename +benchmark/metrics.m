function m = metrics(prediction)
%METRICS Participant-level denominators, abstention and within-ET association.
m=struct('subjects',numel(prediction),'controls',0,'ET',0,'controlFalsePositives',0, ...
    'controlFPR',NaN,'controlFalseSevere',NaN,'clinicalCoverage',NaN, ...
    'controlCoverage',NaN,'ETCoverage',NaN,'mildCoverage',NaN,'severeCoverage',NaN, ...
    'controlIndeterminate',0,'ETIndeterminate',0,'mildIndeterminate',0,'severeIndeterminate',0, ...
    'mildSubjects',0,'mildCaptured',0,'mildCapture',NaN,'moderateCapture',NaN, ...
    'severeSubjects',0,'severeCaptured',0,'severeCapture',NaN, ...
    'severeUndercalled',0,'falseSevereLower',NaN,'ETCapture',NaN, ...
    'balancedLogLoss',NaN,'brier',NaN,'AUCavailable',NaN, ...
    'mildAUCavailable',NaN,'severeAUCavailable',NaN,'severitySpearman',NaN, ...
    'severityKendallTauB',NaN,'severityMAE',NaN,'severityRPS',NaN,'severityPairs',0);
if isempty(prediction), return; end
c=strcmp({prediction.group},'control'); e=strcmp({prediction.group},'ET');
ok=strcmp({prediction.status},'ok'); d=[prediction.decision]; p=[prediction.score];
g=[prediction.trueSeverityGroup]; mild=e&g==1; severe=e&g==3;
m.controls=sum(c); m.ET=sum(e); m.clinicalCoverage=mean(ok);
m.controlCoverage=divide(sum(c&ok),sum(c)); m.ETCoverage=divide(sum(e&ok),sum(e));
m.mildCoverage=divide(sum(mild&ok),sum(mild)); m.severeCoverage=divide(sum(severe&ok),sum(severe));
m.controlIndeterminate=sum(c&~ok); m.ETIndeterminate=sum(e&~ok);
m.mildIndeterminate=sum(mild&~ok); m.severeIndeterminate=sum(severe&~ok);
m.controlFalsePositives=sum(c&d==1); m.controlFPR=divide(m.controlFalsePositives,sum(c));
m.controlFalseSevere=divide(sum(c&[prediction.severityDecision]==3),sum(c));
m.ETCapture=divide(sum(e&d==1),sum(e));
m.mildSubjects=sum(mild); m.mildCaptured=sum(mild&d==1); m.mildCapture=divide(m.mildCaptured,m.mildSubjects);
m.moderateCapture=divide(sum(e&g==2&d==1),sum(e&g==2));
m.severeSubjects=sum(severe); m.severeCaptured=sum(severe&d==1); m.severeCapture=divide(m.severeCaptured,m.severeSubjects);
m.severeUndercalled=sum(severe&ok&[prediction.severityDecision]<3);
m.falseSevereLower=divide(sum(e&g>0&g<3&[prediction.severityDecision]==3),sum(e&g>0&g<3));
y=double(e); p=max(1e-12,min(1-1e-12,p)); loss=-(y.*log(p)+(1-y).*log1p(-p));
if any(c)&&any(e), m.balancedLogLoss=.5*(mean(loss(c))+mean(loss(e))); end
m.brier=mean((p-y).^2);
m.AUCavailable=auc(p(e&ok),p(c&ok));
m.mildAUCavailable=auc(p(mild&ok),p(c&ok)); m.severeAUCavailable=auc(p(severe&ok),p(c&ok));
S=[prediction.trueSeverityScore]; estimate=[prediction.severityPrediction];
use=e&ok&isfinite(S)&isfinite(estimate); m.severityPairs=sum(use);
if sum(use)>=3
    a=tiedRank(S(use)); b=tiedRank(estimate(use)); m.severitySpearman=correlation(a,b);
    m.severityKendallTauB=tau(S(use),estimate(use)); m.severityMAE=mean(abs(S(use)-estimate(use)));
end
use=e&ok&g>0;
if any(use)
    P1=[prediction(use).ordinalP1]; P2=P1+[prediction(use).ordinalP2]; truth=g(use);
    m.severityRPS=mean(((P1-double(truth<=1)).^2+(P2-double(truth<=2)).^2)/2);
end
end
function v=divide(a,b)
v=NaN; if b>0, v=a/b; end
end
function v=auc(positive,negative)
v=NaN; if isempty(positive)||isempty(negative), return; end
D=bsxfun(@minus,positive(:),negative(:)'); v=mean(double(D(:)>0)+.5*double(D(:)==0));
end
function r=tiedRank(x)
[v,ix]=sort(x); r=zeros(size(x)); k=1;
while k<=numel(v)
    j=k; while j<numel(v)&&v(j+1)==v(k), j=j+1; end
    r(ix(k:j))=(k+j)/2; k=j+1;
end
end
function v=correlation(a,b)
a=a-mean(a); b=b-mean(b); v=NaN; q=sqrt(sum(a.^2)*sum(b.^2)); if q>0, v=sum(a.*b)/q; end
end
function v=tau(a,b)
concord=0; discord=0; tieA=0; tieB=0;
for i=1:numel(a)-1
    for j=i+1:numel(a)
        x=sign(a(i)-a(j)); y=sign(b(i)-b(j));
        if x*y>0, concord=concord+1; elseif x*y<0, discord=discord+1;
        elseif x==0&&y~=0, tieA=tieA+1; elseif y==0&&x~=0, tieB=tieB+1; end
    end
end
v=NaN; denominator=sqrt((concord+discord+tieA)*(concord+discord+tieB));
if denominator>0, v=(concord-discord)/denominator; end
end
