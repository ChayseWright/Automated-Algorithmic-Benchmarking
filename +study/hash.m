function value = hash(input)
%HASH SHA-256 of canonical JSON metadata or raw uint8 bytes.
if ~isa(input,'uint8'), input=uint8(unicode2native(jsonencode(canonical(input)),'UTF-8')); end
if exist('OCTAVE_VERSION','builtin')~=0
    value=lower(builtin('hash','sha256',char(input(:).'))); return;
end
if ~usejava('jvm'), error('tbench:Hash','SHA-256 cache validation requires the MATLAB JVM.'); end
engine=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
if ~isempty(input), engine.update(typecast(input(:),'int8')); end
value=hexDigest(engine.digest());
end
function out=canonical(v)
if isa(v,'function_handle'), out=func2str(v);
elseif isnumeric(v)
    if ~isreal(v)
        out=struct('numericClass',class(v),'shape',size(v),'real',canonical(real(v)),'imag',canonical(imag(v)));
    elseif any(~isfinite(v(:)))
        entries=num2cell(v(:)');
        for k=1:numel(entries)
            if isnan(entries{k}), entries{k}='NaN';
            elseif entries{k}==Inf, entries{k}='+Inf';
            elseif entries{k}==-Inf, entries{k}='-Inf'; end
        end
        out=struct('numericClass',class(v),'shape',size(v),'entries',{entries});
    else
        out=struct('numericClass',class(v),'shape',size(v),'values',v);
    end
elseif isstruct(v)
    out=orderfields(v);
    names=fieldnames(out);
    for k=1:numel(out)
        for j=1:numel(names), out(k).(names{j})=canonical(out(k).(names{j})); end
    end
elseif iscell(v)
    out=cell(size(v)); for k=1:numel(v), out{k}=canonical(v{k}); end
else, out=v;
end
end
function value=hexDigest(bytes)
value=lower(reshape(dec2hex(typecast(bytes,'uint8'),2).',1,[]));
end
