function [groups,cuts,description] = severityGroups(subjects,cfg,training)
%SEVERITYGROUPS Never infer clinical grades from an unspecified percentile.
if nargin<3, training=true(1,numel(subjects)); end
score=[subjects.severity_score]; et=strcmp({subjects.group},'ET');
groups=zeros(1,numel(subjects)); cuts=[];
if strcmp(cfg.clinical.severityMode,'provided')
    names={'mild','moderate','severe'};
    for k=1:numel(names), groups(et&strcmp({subjects.severity_group},names{k}))=k; end
    description='provided_clinical_categories';
else
    cuts=cfg.clinical.severityCuts;
    if isempty(cuts)
        v=score(et&training&isfinite(score));
        if isempty(v), cuts=[]; description='severity_unavailable'; return; end
        cuts=tbench.quantile(v,[1/3 2/3]);
        description='training_relative_tertiles';
    else
        if numel(cuts)~=2||any(~isfinite(cuts))||cuts(1)>cuts(2), error('tbench:SeverityCuts','Provide two ordered finite cuts.'); end
        description='prespecified_score_cuts';
    end
    use=et&isfinite(score);
    groups(use&score<=cuts(1))=1;
    groups(use&score>cuts(1)&score<=cuts(2))=2;
    groups(use&score>cuts(2))=3;
end
end
