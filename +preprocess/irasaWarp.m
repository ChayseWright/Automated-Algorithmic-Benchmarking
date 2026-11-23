function out = irasaWarp(y,factor,fs,opts)
%IRASAWARP Anti-aliased down-warp or guarded sinc up-warp for IRASA only.
% The caller intentionally interprets the warped samples at the ORIGINAL Fs.
% This shifts spectral structure; it is not the reporting-domain rate reducer.
if ~isscalar(factor)||~isfinite(factor)||factor<=0
    error('tbench:IRASAFactor','The IRASA warp factor must be finite and positive.');
end
y=y(:);
if factor<=1
    target=fs*factor;
    [out,~]=tbench.preprocess.rateChange(y,fs,target,.40*target,.48*target,opts);
    return;
end
position=(0:floor((numel(y)-1)*factor))'/factor+1;
L=opts.sincHalfLength; center=floor(position); indices=bsxfun(@plus,center,-L+1:L);
distance=bsxfun(@minus,position,indices); s=ones(size(distance)); q=distance~=0;
s(q)=sin(pi*distance(q))./(pi*distance(q));
w=s.*besseli(0,8.6*sqrt(max(0,1-(distance/L).^2)))/besseli(0,8.6);
w(abs(distance)>L)=0; w=bsxfun(@rdivide,w,sum(w,2));
inside=all(indices>=1&indices<=numel(y),2); indices=max(1,min(numel(y),indices));
v=y(indices); good=inside&all(isfinite(v)|abs(w)<1e-12,2); v(~isfinite(v))=0;
value=sum(v.*w,2); out=nan(numel(position),1); out(good)=value(good);
end
