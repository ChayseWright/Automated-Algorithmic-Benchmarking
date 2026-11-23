function result = collect(data,candidates,subjects,indices,selected,cfg,outer,training,waveRidge)
%COLLECT Only outer-test subjects, with their outer-training selected candidate.
result=struct('agreement',struct([]),'waveformPrediction',struct([]));
[grades,~,description]=tbench.benchmark.severityGroups(subjects,cfg,training);
for s=indices(:)'
    k=selected; if numel(selected)>1, k=selected(s); end
    file=load(data.recordFiles{s,k},'payload'); if ~isfield(file.payload,'motion'), continue; end
    for name={'agreement','waveformPrediction'}
        key=name{1}; rows=file.payload.motion.(key);
        if strcmp(key,'waveformPrediction')&&nargin>=9&&~isempty(rows), rows=rows([rows.ridge]==waveRidge); end
        for j=1:numel(rows)
            rows(j).group=subjects(s).group; rows(j).phenotype=subjects(s).phenotype;
            rows(j).severityScore=subjects(s).severity_score; rows(j).severityGroup=grades(s);
            rows(j).severityDefinition=description; rows(j).outerFold=outer;
        end
        result.(key)=append(result.(key),rows);
    end
end
end
function a=append(a,b)
if isempty(b), return; end
if isempty(a), a=b(:); else, a=[a(:);b(:)]; end
end
