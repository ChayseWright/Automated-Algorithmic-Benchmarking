function rows = coherence(packet,ref,target,m)
%COHERENCE Vector-output association and conditional contribution accounting.
% Surrogates permute whole windows, preserving taper bundles and local spectra.
blank=struct('target',target.name,'unit',target.unit,'muscle',0,'status','unavailable', ...
    'frequency',NaN,'coherentFraction',NaN,'peakFraction',NaN,'peakFrequency',NaN, ...
    'multipleFraction',NaN,'isolatedContribution',NaN,'excludedContribution',NaN, ...
    'directPower',NaN,'interferencePower',NaN,'modelSourcePower',NaN,'measuredMotionPower',NaN, ...
    'inputRank',0,'conditionNumber',Inf,'segments',0,'nominalReplicates',0, ...
    'linearModelConditionsMet',false,'attributionStatus','unavailable', ...
    'powerScope','one_sided_power_spectral_density_at_frequency', ...
    'maximumSearchP',NaN,'nullMethod','not_requested','surrogates',m.coherence.surrogates);
rows=blank;
if numel(target.muscles)>1&&~tbench.option(m,'sourceSynchronized',false)
    rows.status='EMG_source_clock_unverified'; return;
end
if isempty(packet)||~all(packet.hasWaveform(target.muscles)), rows.status='native_waveform_unavailable'; return; end
if ~ref.clock.phaseReliable, rows.status='phase_alignment_unavailable'; return; end
M=target.muscles; X=packet.waveform(:,M); X(~packet.valid(:,M))=NaN;
[X,~]=tbench.motion.align(packet.time,X,target,m); Y=target.waveform;
valid=target.valid&all(isfinite(X),2)&all(isfinite(Y),2);
fs=target.fs; N=round(m.coherence.windowSeconds*fs); W=tbench.tapers(N,'dpss',m.coherence.NW,m.coherence.nTapers);
K=size(W,2); nfft=2^nextpow2(N); f=(0:nfft/2)'*fs/nfft;
band=find(f>=m.searchBand(1)&f<=m.searchBand(2)); zX=[]; zY=[]; segments=0;
for a=1:N:size(X,1)-N+1
    ix=a:a+N-1; if ~all(valid(ix)), continue; end
    segments=segments+1;
    for taper=1:K
        vx=bsxfun(@minus,X(ix,:),mean(X(ix,:),1)); vy=bsxfun(@minus,Y(ix,:),mean(Y(ix,:),1));
        x=fft(bsxfun(@times,vx,W(:,taper)),nfft,1)/sqrt(fs);
        y=fft(bsxfun(@times,vy,W(:,taper)),nfft,1)/sqrt(fs);
        oneSided=sqrt(2)*ones(numel(band),1); oneSided(f(band)==0|f(band)==fs/2)=1;
        zX(:,:,segments,taper)=bsxfun(@times,x(band,:),oneSided);
        zY(:,:,segments,taper)=bsxfun(@times,y(band,:),oneSided); %#ok<AGROW>
    end
end
blank.segments=segments; blank.nominalReplicates=segments*K;
if segments<m.coherence.minimumSegments, blank.status='insufficient_independent_windows'; rows=blank; return; end
ordinary=fraction(zX,zY); null=[];
if m.coherence.surrogates>0
    state=rng; cleanup=onCleanup(@()rng(state)); %#ok<NASGU>
    rng(m.coherence.seed+segments*17+numel(M),'twister'); null=zeros(m.coherence.surrogates,1);
    for b=1:numel(null)
        order=randperm(segments); value=fraction(zX,zY(:,:,order,:)); finite=value(isfinite(value));
        if isempty(finite), null(b)=NaN; else, null(b)=max(finite); end
    end
    blank.nullMethod='whole_window_pairing_max_over_frequencies_and_inputs_exchangeability_required';
end
motionPower=zeros(numel(band),1);
for j=1:numel(band), y=reshape(zY(j,:,:,:),size(Y,2),[])'; motionPower(j)=mean(sum(abs(y).^2,2)); end
if isfinite(target.frequency), [~,at]=min(abs(f(band)-target.frequency));
else, [~,at]=max(motionPower);
end
A=reshape(zX(at,:,:,:),numel(M),[])'; B=reshape(zY(at,:,:,:),size(Y,2),[])';
Sxx=(A.'*conj(A))/size(A,1); Sxy=(A.'*conj(B))/size(A,1); Syy=(B.'*conj(B))/size(B,1);
Sxx=(Sxx+Sxx')/2; Syx=Sxy'; output=real(trace(Syy));
e=real(eig(Sxx)); tolerance=max(max(e),realmin)*1e-8; inputRank=sum(e>tolerance);
condition=Inf; if min(e)>0, condition=max(e)/min(e); end
full=inputRank==numel(M)&&condition<=m.coherence.maximumCondition&&segments*K>=m.coherence.replicatesPerInput*numel(M);
if full, H=Syx/Sxx;
else, H=Syx/(Sxx+eye(numel(M))*max(real(trace(Sxx))/numel(M),realmin)*m.coherence.shrinkage);
end
modeled=real(trace(H*Sxx*H')); multiple=NaN;
if output>0, multiple=real(trace(H*Sxy))/output; end
eligible=full&&isfinite(multiple)&&multiple>=m.coherence.minimumMultipleCoherence;
direct=zeros(numel(M),1); interference=zeros(numel(M));
for a=1:numel(M)
    direct(a)=real(Sxx(a,a))*sum(abs(H(:,a)).^2);
    for b=a+1:numel(M)
        value=2*real(Sxx(a,b)*sum(H(:,a).*conj(H(:,b)))); interference(a,b)=value; interference(b,a)=value;
    end
end
rows=repmat(blank,1,numel(M));
for j=1:numel(M)
    r=blank; r.muscle=M(j); r.status='ok'; r.frequency=f(band(at)); r.coherentFraction=ordinary(at,j);
    v=ordinary(:,j); finite=find(isfinite(v));
    if ~isempty(finite), [r.peakFraction,pick]=max(v(finite)); r.peakFrequency=f(band(finite(pick))); end
    r.multipleFraction=multiple; r.inputRank=inputRank; r.conditionNumber=condition;
    r.modelSourcePower=modeled; r.measuredMotionPower=output;
    r.linearModelConditionsMet=eligible; r.attributionStatus='insufficient_rank_replicates_or_model_fit';
    if eligible
        r.directPower=direct(j); r.interferencePower=sum(interference(j,:));
        r.isolatedContribution=direct(j)/output;
        r.excludedContribution=(direct(j)+sum(interference(j,:)))/output;
        r.attributionStatus='conditional_linear_model_not_interventional_causality';
    end
    good=null(isfinite(null));
    if ~isempty(good)&&isfinite(r.peakFraction), r.maximumSearchP=(1+sum(good>=r.peakFraction))/(1+numel(good)); end
    rows(j)=r;
end
end
function out=fraction(zX,zY)
F=size(zX,1); M=size(zX,2); D=size(zY,2); out=nan(F,M);
for j=1:F
    X=reshape(zX(j,:,:,:),M,[])'; Y=reshape(zY(j,:,:,:),D,[])';
    cross=X.'*conj(Y)/size(X,1); px=mean(abs(X).^2,1)'; py=mean(sum(abs(Y).^2,2));
    good=px>0&py>0; out(j,good)=sum(abs(cross(good,:)).^2,2)./(px(good)*py);
end
end
