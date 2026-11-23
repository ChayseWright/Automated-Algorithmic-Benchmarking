function model = fitLinear(X,y,lambda,weights)
%FITLINEAR Small signed linear readout; scaling and variance are training-only.
if nargin<4||isempty(weights), weights=ones(size(y)); end
good=all(isfinite(X),2)&isfinite(y)&isfinite(weights)&weights>0;
X=X(good,:); y=y(good); weights=weights(good);
if isempty(y), model=[]; return; end
weights=weights/mean(weights); mu=mean(X,1); sd=std(X,1,1); sd(sd<1e-10)=1;
A=[ones(size(X,1),1) bsxfun(@rdivide,bsxfun(@minus,X,mu),sd)];
penalty=lambda*eye(size(A,2)); penalty(1,1)=0;
beta=(A'*bsxfun(@times,A,weights)+penalty)\(A'*(weights.*y)); residual=y-A*beta;
variance=max(1e-8,sum(weights.*residual.^2)/sum(weights));
model=struct('mu',mu,'sd',sd,'beta',beta,'variance',variance,'trainingRows',numel(y),'lambda',lambda);
end
