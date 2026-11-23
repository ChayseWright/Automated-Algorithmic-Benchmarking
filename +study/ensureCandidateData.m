function data = ensureCandidateData(data,subjects,subjectIndices,candidateIndices,candidates,cfg,codeHash)
%ENSURECANDIDATEDATA Fill requested subjects only; physical caches are reusable.
needed=[];
for k=candidateIndices(:)'
    for s=subjectIndices(:)'
        if k>size(data.recordFiles,2)||isempty(data.recordFiles{s,k})
            needed(end+1)=k; break; %#ok<AGROW>
        end
    end
end
needed=unique(needed,'stable'); if isempty(needed), return; end
addition=tbench.study.extract(subjects(subjectIndices),candidates(needed),cfg,codeHash);
data=tbench.study.mergeCandidateData(data,addition,subjectIndices,needed,cfg);
end
