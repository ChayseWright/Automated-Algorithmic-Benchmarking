function rows = qualitySummary(F,cfg)
%QUALITYSUMMARY B8 all enrolled workflow counts, by task and clinical group.
rows=struct([]); if isempty(F), return; end
for group={'control','ET'}
    for j=1:numel(cfg.tasks)
        Q=F(strcmp({F.group},group{1})&strcmp({F.task},cfg.tasks{j}));
        if isempty(Q), continue; end
        statuses=unique({Q.status});
        for k=1:numel(statuses)
            J=Q(strcmp({Q.status},statuses{k}));
            r=struct('pipelineId',J(1).pipelineId,'group',group{1},'task',cfg.tasks{j}, ...
                'status',statuses{k},'records',numel(J),'participants',numel(unique({J.subject})), ...
                'meanValidFraction',mean([J.validFraction]),'meanReviewFraction',mean([J.reviewFraction]), ...
                'broadWitnessOutsidePrimary',sum([J.guardOutsidePrimary]));
            if isempty(rows), rows=r; else, rows(end+1)=r; end %#ok<AGROW>
        end
    end
end
end
