function value = sourceHash(spec)
%SOURCEHASH Retain identity of raw content, chosen root and MVC reference.
if isfield(spec,'mvc')&&~isempty(spec.mvc)
    v=spec.mvc;
    if ~isnumeric(v)||~isreal(v)||numel(v)~=15||any(~isnan(v(:))&(~isfinite(v(:))|v(:)<=0))
        error('tbench:MVC','MVC overrides require 15 positive finite amplitudes or NaNs.');
    end
end
if isfield(spec,'file')&&~isempty(spec.file)
    raw=tbench.study.fileHash(spec.file);
else
    loaded=tbench.loadSubject(spec); tasks=sort(fieldnames(loaded.data)); parts=cell(1,0);
    tasks=tasks(~cellfun('isempty',regexp(tasks,'^[PK][1-7]$','once')));
    for k=1:numel(tasks)
        if ~isstruct(loaded.data.(tasks{k}))
            parts{end+1}=struct('task',tasks{k},'trial','invalid_task','size',size(loaded.data.(tasks{k})), ...
                'sha',tbench.study.hash(loaded.data.(tasks{k}))); %#ok<AGROW>
            continue;
        end
        trials=sort(fieldnames(loaded.data.(tasks{k})));
        for j=1:numel(trials)
            A=loaded.data.(tasks{k}).(trials{j});
            isComplex=false;
            if isnumeric(A)
                values=full(double(A(:)));
                isComplex=~isreal(values);
                if isComplex, values=[real(values);imag(values)]; end
                digest=tbench.study.hash(typecast(values,'uint8'));
            else, digest=tbench.study.hash(A);
            end
            parts{end+1}=struct('task',tasks{k},'trial',trials{j},'size',size(A), ...
                'sha',tbench.study.hash(struct('complex',isComplex,'digest',digest))); %#ok<AGROW>
        end
    end
    raw=tbench.study.hash(parts);
end
variable=tbench.option(spec,'variable','');
mvc=tbench.option(spec,'mvc',[]);
value=tbench.study.hash(struct('raw',raw,'variable',variable,'mvc',mvc));
end
