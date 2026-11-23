function model = fitLogistic(X,y,lambda,balanced)
%FITLOGISTIC Small ridge model, training-only imputation/scale, stable IRLS.
if nargin<4, balanced=true; end
y=double(y(:)); n=numel(y); p=size(X,2);
if n==0, error('tbench:Readout','No training observations.'); end
fill=zeros(1,p);
for k=1:p
    v=X(:,k); v=v(isfinite(v)); if ~isempty(v), fill(k)=median(v); end
    X(~isfinite(X(:,k)),k)=fill(k);
end
mu=mean(X,1); sd=std(X,1,1); sd(sd<1e-10)=1;
A=[ones(n,1) bsxfun(@rdivide,bsxfun(@minus,X,mu),sd)];
prior=(sum(y)+.5)/(n+1); b=zeros(p+1,1); b(1)=log(prior/(1-prior));
w=ones(n,1);
if balanced&&any(y==0)&&any(y==1)
    w(y==0)=n/(2*sum(y==0)); w(y==1)=n/(2*sum(y==1));
end
penalty=lambda*eye(p+1); penalty(1,1)=0;
constant=~any(y==0)||~any(y==1);
if ~constant
    for it=1:100
        eta=A*b; probability=1./(1+exp(-max(-40,min(40,eta))));
        H=A'*(bsxfun(@times,A,w.*max(probability.*(1-probability),1e-8)))+penalty+eye(p+1)*1e-10;
        gradient=A'*(w.*(probability-y))+penalty*b; step=H\gradient;
        before=objective(A*b,y,w,b,penalty); rate=1;
        while rate>1e-6
            next=b-rate*step;
            if objective(A*next,y,w,next,penalty)<=before+1e-12, break; end
            rate=rate/2;
        end
        b=next; if norm(rate*step)<1e-7, break; end
    end
end
model=struct('fill',fill,'mu',mu,'sd',sd,'beta',b,'prior',prior, ...
    'lambda',lambda,'constant',constant,'trainingRows',n);
end
function v=objective(eta,y,w,b,penalty)
soft=max(eta,0)+log1p(exp(-abs(eta))); v=sum(w.*(soft-y.*eta))+.5*b'*penalty*b;
end
