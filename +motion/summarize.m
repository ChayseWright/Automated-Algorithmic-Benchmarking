function [summary,retention,contrasts] = summarize(predictions,cfg)
%SUMMARIZE Operational mechanical retention with paired subject uncertainty.
summary=struct([]); retention=struct([]); contrasts=struct([]);
if isempty(predictions), return; end
state=rng; cleanup=onCleanup(@()rng(state)); rng(cfg.seed+8401,'twister'); %#ok<NASGU>
ids=unique({predictions.pipelineId}); levels=unique({predictions.level});
for a=1:numel(ids)
    for k=1:numel(levels)
        Q=predictions(strcmp({predictions.pipelineId},ids{a})&strcmp({predictions.level},levels{k}));
        if isempty(Q), continue; end
        subsets={'all','control','ET','lower_severity','higher_severity'};
        for j=1:numel(subsets)
            R=subset(Q,subsets{j}); if isempty(R), continue; end
            s=tbench.motion.utility(R,cfg); s.pipelineId=ids{a}; s.level=levels{k}; s.subgroup=subsets{j};
            summary=append(summary,s);
        end
    end
    rich=predictions(strcmp({predictions.pipelineId},ids{a})&strcmp({predictions.level},'R'));
    for k=1:numel(levels)
        if strcmp(levels{k},'R'), continue; end
        Q=predictions(strcmp({predictions.pipelineId},ids{a})&strcmp({predictions.level},levels{k}));
        if isempty(Q)||isempty(rich), continue; end
        for subgroup={'all','control','ET','lower_severity','higher_severity'}
            [C,R]=align(subset(Q,subgroup{1}),subset(rich,subgroup{1})); if isempty(C), continue; end
            q=tbench.motion.utility(C,cfg); r=tbench.motion.utility(R,cfg); retained=NaN;
            if r.baselineLoss-r.loss>cfg.motion.minimumBaselineLoss, retained=(r.baselineLoss-q.loss)/(r.baselineLoss-r.loss); end
            [ci,estimable]=retentionCI(C,R,cfg);
            row=struct('pipelineId',ids{a},'level',levels{k},'subgroup',subgroup{1},'subjects',q.subjects, ...
                'richLoss',r.loss,'compressedLoss',q.loss,'baselineLoss',r.baselineLoss, ...
                'lossIncrease',q.loss-r.loss,'retainedUtility',retained,'CI025',ci(1),'CI975',ci(2), ...
                'bootstrapEstimable',estimable,'bootstrapSamples',cfg.tuning.bootstrapSamples, ...
                'uncertaintyScope','subject_bootstrap_conditional_on_outer_fits', ...
                'scope','held_out_RMS_amplitude_decoder_not_total_information');
            retention=append(retention,row);
        end
    end
end
% Paired subject bootstrap of already fitted outer models, not full re-nesting.
for a=1:numel(ids)-1
    for b=a+1:numel(ids)
        A=predictions(strcmp({predictions.pipelineId},ids{a})&strcmp({predictions.level},'F1'));
        B=predictions(strcmp({predictions.pipelineId},ids{b})&strcmp({predictions.level},'F1'));
        [A,B]=align(A,B); if isempty(A), continue; end
        subjects=unique({A.subject}); differences=nan(numel(subjects),1);
        for k=1:numel(subjects)
            sa=tbench.motion.utility(A(strcmp({A.subject},subjects{k})),cfg);
            sb=tbench.motion.utility(B(strcmp({B.subject},subjects{k})),cfg); differences(k)=sa.gain-sb.gain;
        end
        differences=differences(isfinite(differences)); if isempty(differences), continue; end
        samples=nan(cfg.tuning.bootstrapSamples,1);
        for z=1:numel(samples), samples(z)=mean(differences(randi(numel(differences),numel(differences),1))); end
        ci=tbench.quantile(samples,[.025 .975]);
        row=struct('pipelineA',ids{a},'pipelineB',ids{b},'subjects',numel(differences), ...
            'gainDifference',mean(differences),'CI025',ci(1),'CI975',ci(2), ...
            'scope','paired_subject_bootstrap_conditional_on_outer_fits');
        contrasts=append(contrasts,row);
    end
end
end
function [ci,estimable]=retentionCI(C,R,cfg)
ci=[NaN NaN]; estimable=0; ids=unique({C.subject}); values=[];
for k=1:numel(ids)
    c=tbench.motion.utility(C(strcmp({C.subject},ids{k})),cfg);
    r=tbench.motion.utility(R(strcmp({R.subject},ids{k})),cfg);
    v=[c.loss r.loss r.baselineLoss];
    if all(isfinite(v)), values(end+1,:)=v; end %#ok<AGROW>
end
if size(values,1)<2||cfg.tuning.bootstrapSamples<1, return; end
ix=randi(size(values,1),size(values,1),cfg.tuning.bootstrapSamples);
c=mean(reshape(values(ix,1),size(ix)),1); r=mean(reshape(values(ix,2),size(ix)),1);
b=mean(reshape(values(ix,3),size(ix)),1); good=b-r>cfg.motion.minimumBaselineLoss;
estimable=sum(good); ratios=(b(good)-c(good))./(b(good)-r(good));
if estimable>=.8*cfg.tuning.bootstrapSamples, ci=tbench.quantile(ratios,[.025 .975]); end
end
function rows=subset(rows,name)
switch name
    case 'control', rows=rows(strcmp({rows.group},'control'));
    case 'ET', rows=rows(strcmp({rows.group},'ET'));
    case 'lower_severity', rows=rows(strcmp({rows.group},'ET')&[rows.severityGroup]==1);
    case 'higher_severity', rows=rows(strcmp({rows.group},'ET')&[rows.severityGroup]==3);
end
end
function [A,B]=align(A,B)
if isempty(A)||isempty(B), A=[]; B=[]; return; end
a=arrayfun(@(r)[r.subject '|' r.task '|' r.trial '|' r.target],A,'UniformOutput',false);
b=arrayfun(@(r)[r.subject '|' r.task '|' r.trial '|' r.target],B,'UniformOutput',false);
[~,ia,ib]=intersect(a,b,'stable'); A=A(ia); B=B(ib);
end
function a=append(a,b)
if isempty(b), return; end
if isempty(a), a=b(:); else, a=[a(:);b(:)]; end
end
