function y = predictLinear(model,X)
%PREDICTLINEAR Unsupported predictors remain NaN for an explicit baseline fallback.
y=nan(size(X,1),1); if isempty(model), return; end
good=all(isfinite(X),2); A=[ones(sum(good),1) bsxfun(@rdivide,bsxfun(@minus,X(good,:),model.mu),model.sd)];
y(good)=A*model.beta;
end
