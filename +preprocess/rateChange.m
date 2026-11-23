function [y,t,audit] = rateChange(x,fs,targetFs,passHz,stopHz,opts)
%RATECHANGE Guarded anti-alias FIR and integer or windowed-sinc resampling.
% Invalid support remains invalid. Each column uses exactly the same clock.
if targetFs>fs*(1+1e-10), error('tbench:Upsampling','Reduction module never invents a higher acquisition rate.'); end
if stopHz>=targetFs/2
    if abs(stopHz-targetFs/2)<1e-9, stopHz=stopHz*(1-1e-6);
    else, error('tbench:OutputNyquist','Stopband must begin no higher than output Nyquist.'); end
end
[h,audit]=tbench.preprocess.designAA(fs,passHz,stopHz,opts.stopAttenuationDB,opts.passRippleDB);
filtered=tbench.applyFIR(x,h); t=(0:floor((size(x,1)-1)/fs*targetFs))'/targetFs;
position=t*fs+1;
if all(abs(position-round(position))<1e-7)
    y=filtered(round(position),:); method='integer_decimation';
    if abs(targetFs-fs)<1e-9, method='same_rate_filter'; end
else
    L=opts.sincHalfLength; center=floor(position); offsets=-L+1:L;
    indices=bsxfun(@plus,center,offsets); distance=bsxfun(@minus,position,indices);
    s=ones(size(distance)); q=distance~=0; s(q)=sin(pi*distance(q))./(pi*distance(q));
    w=besseli(0,8.6*sqrt(max(0,1-(distance/L).^2)))/besseli(0,8.6);
    w=w.*s; w(abs(distance)>L)=0; w=bsxfun(@rdivide,w,sum(w,2));
    inside=all(indices>=1&indices<=size(x,1),2); indices=max(1,min(size(x,1),indices));
    y=nan(numel(t),size(x,2));
    for m=1:size(x,2)
        v=filtered(indices+(m-1)*size(x,1)); valid=inside&all(isfinite(v)|abs(w)<1e-12,2);
        v(~isfinite(v))=0; value=sum(v.*w,2); y(valid,m)=value(valid);
    end
    audit.halfSupportSeconds=audit.halfSupportSeconds+L/fs; method='windowed_sinc';
end
audit.outputFs=targetFs; audit.ratio=fs/targetFs; audit.method=method;
audit.interpolationHalfLength=opts.sincHalfLength;
end
