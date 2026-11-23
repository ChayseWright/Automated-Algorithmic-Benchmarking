function P = frozenCompression(rich,compressed,subjects,train,test,cfg,model,threshold)
%FROZENCOMPRESSION Remove descriptors without refitting the rich clinical model.
% Omitted raw descriptors are fixed at TRAINING medians; mask differences abstain.
reference=rich.X; fill=zeros(1,size(reference,2)); train=train(:);
for j=1:size(reference,2)
    v=reference(train&rich.available,j); v=v(isfinite(v));
    if ~isempty(v), fill(j)=median(v); end
end
X=repmat(fill,size(reference,1),1);
[present,position]=ismember(compressed.featureNames,rich.featureNames);
if ~all(present), error('tbench:CompressionFeatures','Compression descriptors must be included in the rich reference.'); end
X(:,position)=compressed.X;
available=rich.available&compressed.available;
P=tbench.benchmark.predictClinical(model,X(test,:),subjects(test),available(test),cfg,threshold);
end
