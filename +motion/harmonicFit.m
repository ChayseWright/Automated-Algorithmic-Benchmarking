function a = harmonicFit(y,t,f,fs)
%HARMONICFIT Signed-coordinate fit and conditional HAC power uncertainty.
% The interval does not correct frequency-search selection bias.
u=t(:)-mean(t); H=min(3,floor(.45*fs/f)); X=[ones(size(u)) u];
for h=1:H, X=[X cos(2*pi*h*f*u) sin(2*pi*h*f*u)]; end %#ok<AGROW>
inverse=pinv(X); B=inverse*y; residual=y-X*B; G=X(:,3:4); D=size(y,2); n=numel(t);
influence=zeros(n,2*D);
for d=1:D, influence(:,2*d-1:2*d)=bsxfun(@times,inverse(3:4,:)',residual(:,d)); end
C=influence'*influence; L=min(floor(n/4),ceil(fs/3));
for lag=1:L
    term=influence(1+lag:end,:)'*influence(1:end-lag,:);
    C=C+(1-lag/(L+1))*(term+term');
end
[U,E]=eig((C+C')/2); C=U*diag(max(0,real(diag(E))))*U';
A=kron(eye(D),(G'*G)/n); beta=reshape(B(3:4,:),[],1);
raw=real(beta'*A*beta); bias=real(trace(A*C));
variance=max(0,4*real(beta'*A*C*A*beta)+2*real(trace(A*C*A*C)));
fundamental=G*B(3:4,:);
a=struct('powerRaw',raw,'powerSigned',raw-bias,'powerSD',sqrt(variance), ...
    'noiseBias',bias,'fundamental',fundamental,'frequency',f, ...
    'intervalScope','conditional_selected_frequency_HAC_not_whole_search', ...
    'HAC_lags',L,'harmonics',H);
end
