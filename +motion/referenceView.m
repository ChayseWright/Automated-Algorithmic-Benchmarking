function view = referenceView(runtime,subjects,cfg)
%REFERENCEVIEW Motion-only clinical comparator; no EMG predictors or MVC use.
N=numel(subjects); X=nan(N,3); available=false(N,1); summary=struct([]);
for s=1:N
    values=nan(numel(cfg.tasks),3); coverage=zeros(numel(cfg.tasks),1);
    for j=1:numel(cfg.tasks)
        p=[]; u=[]; d=[]; c=[];
        for r=1:numel(cfg.trials)
            file=runtime.files{s,j,r}; loaded=load(file,'payload'); ref=loaded.payload;
            for k=1:numel(ref.targets)
                t=ref.targets(k);
                if ~isempty(cfg.motion.primaryTargets)&&~any(strcmp(cfg.motion.primaryTargets,t.name)), continue; end
                % One coordinate kind/unit only; never add mm^2 and rad^2.
                if ~strcmp(t.kind,cfg.motion.clinicalTargetKind), continue; end
                if isfinite(t.meanBandPower), p(end+1)=t.meanBandPower; u(end+1)=t.upperWindowPower; d(end+1)=t.positiveSeconds>0; c(end+1)=t.validFraction; end %#ok<AGROW>
            end
        end
        u=u(isfinite(u));
        if ~isempty(p)&&~isempty(u), values(j,:)=[mean(p) max(u) mean(d)]; coverage(j)=mean(c); end
    end
    post=cellfun(@(x)x(1)=='P',cfg.tasks); kin=~post;
    vp=values(post,:); vk=values(kin,:);
    if any(all(isfinite(vp),2))&&any(all(isfinite(vk),2))
        X(s,:)=.5*(mean(vp(all(isfinite(vp),2),:),1)+mean(vk(all(isfinite(vk),2),:),1));
        available(s)=mean(coverage)>=cfg.motion.minimumReferenceCoverage;
    end
    row=struct('subject',subjects(s).id,'group',subjects(s).group,'level','motion_F1', ...
        'coverage',mean(coverage),'status','indeterminate','phenotype',subjects(s).phenotype);
    if available(s), row.status='ok'; end
    if isempty(summary), summary=row; else, summary(end+1)=row; end %#ok<AGROW>
end
view=struct('X',X,'available',available,'summary',summary, ...
    'featureNames',{{'meanMotionBandPower','upperMotionPower','rhythmYield'}});
end
