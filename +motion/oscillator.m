function out = oscillator(y,valid,fs,m)
%OSCILLATOR Motion-only switching bank of damped linear harmonic filters.
% Approximate Gaussian switching evidence, not calibrated posterior confidence.
f=(m.searchBand(1):.25:m.searchBand(2))'; F=numel(f); D=size(y,2); N=size(y,1);
out=struct('frequency',nan(N,1),'onScore',nan(N,1),'model','switching_motion_oscillator_bank');
good=valid&all(isfinite(y),2); if sum(good)<20, return; end
% One physical vector scale preserves rotation invariance. Per-axis whitening
% would promote a quiet/noisy axis to the amplitude of the trembling axis.
scale=max(sqrt(sum(var(y(good,:),1,1))/D),1e-12);
z=y/scale; offVar=max(1,mean(var(z(good,:),1,1)));
q=m.oscillator.processFraction; rho=m.oscillator.damping;
c=rho*cos(2*pi*f/fs); s=rho*sin(2*pi*f/fs);
x=zeros(F,D); v=x; p00=ones(F,1); p01=zeros(F,1); p11=p00;
w=[.9;repmat(.1/F,F,1)];
for n=1:N
    if ~good(n)
        x(:)=0; v(:)=0; p00(:)=1; p01(:)=0; p11(:)=1; w=[.9;repmat(.1/F,F,1)]; continue;
    end
    xp=bsxfun(@times,c,x)-bsxfun(@times,s,v); vp=bsxfun(@times,s,x)+bsxfun(@times,c,v);
    a=c.^2.*p00-2*c.*s.*p01+s.^2.*p11+q;
    b=c.*s.*(p00-p11)+(c.^2-s.^2).*p01;
    d=s.^2.*p00+2*c.*s.*p01+c.^2.*p11+q;
    innovation=bsxfun(@minus,z(n,:),xp); variance=a+.15;
    k0=a./variance; k1=b./variance;
    x=xp+bsxfun(@times,k0,innovation); v=vp+bsxfun(@times,k1,innovation);
    p00=(1-k0).*a; p01=(1-k0).*b; p11=max(realmin,d-k1.*b);
    prior=[w(1)*m.oscillator.stayOff+sum(w(2:end))*(1-m.oscillator.stayOn); ...
        w(1)*(1-m.oscillator.stayOff)/F+m.oscillator.stayOn*(.99*w(2:end)+.01*sum(w(2:end))/F)];
    likelihood=-.5*(D*log(variance)+sum(innovation.^2,2)./variance);
    off=-.5*(D*log(offVar)+sum(z(n,:).^2)/offVar);
    logWeight=log(max(prior,realmin))+[off;likelihood]; logWeight=logWeight-max(logWeight);
    w=exp(max(-700,logWeight)); w=w/sum(w);
    out.onScore(n)=1-w(1);
    if sum(w(2:end))>1e-10, out.frequency(n)=sum(f.*w(2:end))/sum(w(2:end)); end
end
end
