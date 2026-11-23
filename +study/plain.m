function out = plain(value)
%PLAIN JSON-compatible provenance; full native configuration also stays in MAT.
if isa(value,'function_handle'), out=func2str(value);
elseif isstruct(value)
    out=value; names=fieldnames(value);
    for k=1:numel(value)
        for j=1:numel(names), out(k).(names{j})=tbench.study.plain(value(k).(names{j})); end
    end
elseif iscell(value)
    out=cell(size(value)); for k=1:numel(value), out{k}=tbench.study.plain(value{k}); end
else, out=value;
end
end
