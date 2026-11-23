function rows = explanatorySummary(detail,kind)
%EXPLANATORYSUMMARY Trial summaries; complete rows remain in streamed CSVs.
rows=struct([]); if isempty(detail), return; end
methods=unique({detail.method});
for a=1:numel(methods)
    F=detail(strcmp({detail.method},methods{a}));
    if strcmp(kind,'prediction')
        levels=unique({F.level});
        for j=1:numel(levels)
            Q=F(strcmp({F.level},levels{j}));
            r=struct('method',methods{a},'level',levels{j}, ...
                'logScoreGain',meanFinite([Q.logScoreGain]),'coverage',meanFinite([Q.coverage]), ...
                'mseSkill',meanFinite([Q.mseSkill]),'targetFoldRows',numel(Q),'testSamples',sum([Q.testSamples]));
            rows=append(rows,r);
        end
    elseif strcmp(kind,'network')
        r=struct('method',methods{a},'referenceMSC',meanFinite([F.referenceMSC]), ...
            'candidateMSC',meanFinite([F.candidateMSC]), ...
            'meanAbsolutePhaseError',meanFinite(abs([F.phaseError])),'pairs',numel(F));
        rows=append(rows,r);
    elseif strcmp(kind,'frequencySupport')
        r=struct('method',methods{a},'nearestReferenceError',meanFinite([F.nearestReferenceError]), ...
            'possibleOctaveFraction',meanFinite([F.possibleOctaveError]),'muscles',numel(F));
        rows=append(rows,r);
    else
        r=struct('method',methods{a},'referenceBandPower',meanFinite([F.referenceBandPower]), ...
            'residualBandPower',meanFinite([F.residualBandPower]), ...
            'absoluteVarianceIdentityError',meanFinite(abs([F.varianceIdentityError])),'muscles',numel(F));
        rows=append(rows,r);
    end
end
end
function v=meanFinite(x)
x=x(isfinite(x)); v=NaN; if ~isempty(x), v=mean(x); end
end
function a=append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
