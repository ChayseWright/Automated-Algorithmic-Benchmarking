function [y,valid] = sample(time,x,query,half)
%SAMPLE Guarded sinc evaluation of an already anti-aliased uniform signal.
% This does not bridge gaps, extrapolate or assert additional bandwidth.
if nargin<4, half=20; end
time=time(:); query=query(:); y=nan(numel(query),size(x,2)); valid=false(size(y));
if numel(time)<2||size(x,1)~=numel(time), return; end
fs=1/median(diff(time)); position=(query-time(1))*fs+1;
if all(abs(position-round(position))<1e-7)
    ix=round(position); inside=ix>=1&ix<=numel(time); y(inside,:)=x(ix(inside),:); valid=isfinite(y); return;
end
indices=bsxfun(@plus,floor(position),-half+1:half); distance=bsxfun(@minus,position,indices);
w=ones(size(distance)); q=distance~=0; w(q)=sin(pi*distance(q))./(pi*distance(q));
w=w.*besseli(0,8.6*sqrt(max(0,1-(distance/half).^2)))/besseli(0,8.6);
w(abs(distance)>half)=0; w=bsxfun(@rdivide,w,sum(w,2));
inside=all(indices>=1&indices<=numel(time),2); indices=max(1,min(numel(time),indices));
for d=1:size(x,2)
    v=x(indices+(d-1)*size(x,1)); good=inside&all(isfinite(v)|abs(w)<1e-12,2);
    v(~isfinite(v))=0; value=sum(v.*w,2); y(good,d)=value(good); valid(:,d)=good;
end
end
