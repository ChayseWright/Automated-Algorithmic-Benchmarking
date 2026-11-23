function rows = rhythmSummary(F,cfg)
%RHYTHMSUMMARY B1 accepted rhythms, separate from a clinical false positive.
rows=struct([]); if isempty(F), return; end
ids=unique({F.subject});
for s=1:numel(ids)
    Q=F(strcmp({F.subject},ids{s})); taskRows=struct([]);
    for j=1:numel(cfg.tasks)
        J=Q(strcmp({Q.task},cfg.tasks{j}));
        good=strcmp({J.status},'ok')&isfinite([J.powerMean]); K=J(good);
        r=blank(Q(1)); r.task=cfg.tasks{j}; r.plannedRecords=numel(cfg.trials)*15;
        r.eligibleRecords=numel(K); r.failedRecords=numel(J)-numel(K);
        r.missingRecords=max(0,r.plannedRecords-numel(J));
        if ~isempty(K)
            detected=[K.detected]; f=[K.frequencyHz];
            r.yieldPrimary=mean(detected&f>=cfg.searchBand(1)&f<=cfg.searchBand(2));
            r.yield4to8=mean(detected&f>=4&f<=8);
            r.meanPower=finiteMean([K.powerMean]); r.meanMVCPower=finiteMean([K.mvcPowerMean]);
            r.occupancy=finiteMean([K.occupancy]); r.activeWindowPower90=finiteMean([K.mvcActiveWindowPower90]);
            r.validMuscleMinutes=sum([K.validSeconds])/60;
            if any([K.dynamic])&&r.validMuscleMinutes>0
                r.episodesPerMuscleMinute=sum([K([K.dynamic]).episodeCount])/r.validMuscleMinutes;
            end
        end
        taskRows=append(taskRows,r); rows=append(rows,r);
    end
    r=blank(Q(1)); r.task='task_balanced';
    r.plannedRecords=sum([taskRows.plannedRecords]); r.eligibleRecords=sum([taskRows.eligibleRecords]);
    r.failedRecords=sum([taskRows.failedRecords]); r.missingRecords=sum([taskRows.missingRecords]);
    r.validMuscleMinutes=sum([taskRows.validMuscleMinutes]);
    post=cellfun(@(x)x(1)=='P',{taskRows.task});
    for name={'yieldPrimary','yield4to8','meanPower','meanMVCPower','occupancy','activeWindowPower90','episodesPerMuscleMinute'}
        field=name{1}; a=finiteMean([taskRows(post).(field)]); b=finiteMean([taskRows(~post).(field)]);
        if isfinite(a)&&isfinite(b), r.(field)=.5*(a+b); end
    end
    rows=append(rows,r);
end
end
function r=blank(F)
r=struct('subject',F.subject,'pipelineId',F.pipelineId,'group',F.group,'phenotype',F.phenotype, ...
    'task','','plannedRecords',0,'eligibleRecords',0,'failedRecords',0,'missingRecords',0, ...
    'yieldPrimary',NaN,'yield4to8',NaN,'meanPower',NaN,'meanMVCPower',NaN, ...
    'occupancy',NaN,'activeWindowPower90',NaN,'episodesPerMuscleMinute',NaN,'validMuscleMinutes',0);
end
function v=finiteMean(x)
x=x(isfinite(x)); v=NaN; if ~isempty(x), v=mean(x); end
end
function a=append(a,b)
if isempty(a), a=b; else, a(end+1)=b; end
end
