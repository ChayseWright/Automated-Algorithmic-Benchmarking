function p = predictLogistic(model,X)
for k=1:size(X,2), X(~isfinite(X(:,k)),k)=model.fill(k); end
A=[ones(size(X,1),1) bsxfun(@rdivide,bsxfun(@minus,X,model.mu),model.sd)];
p=1./(1+exp(-max(-40,min(40,A*model.beta))));
end
