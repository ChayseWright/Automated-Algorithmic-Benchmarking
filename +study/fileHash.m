function value = fileHash(filename)
%FILEHASH Stream file bytes; source content, not its name/mtime, keys the cache.
fid=fopen(filename,'rb'); if fid<0, error('tbench:Source','Cannot read %s.',filename); end
cleanup=onCleanup(@()fclose(fid)); %#ok<NASGU>
if exist('OCTAVE_VERSION','builtin')~=0
    bytes=fread(fid,Inf,'*uint8'); value=tbench.study.hash(bytes); return;
end
if ~usejava('jvm'), error('tbench:Hash','Source hashing requires the MATLAB JVM.'); end
engine=javaMethod('getInstance','java.security.MessageDigest','SHA-256');
while ~feof(fid)
    bytes=fread(fid,1024*1024,'*uint8'); if isempty(bytes), break; end
    engine.update(typecast(bytes(:),'int8'));
end
value=lower(reshape(dec2hex(typecast(engine.digest(),'uint8'),2).',1,[]));
end
