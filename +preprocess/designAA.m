function [h,audit] = designAA(fs,passHz,stopHz,attenuation,rippleDB)
%DESIGNAA Kaiser-windowed sinc with measured pass/stop checks, base MATLAB.
if ~(0<passHz&&passHz<stopHz&&stopHz<fs/2)
    error('tbench:AAEdges','Require 0 < pass < stop < input Nyquist.');
end
if nargin<4, attenuation=90; end
if nargin<5, rippleDB=.05; end
persistent keys coefficients reports
if isempty(keys), keys={}; coefficients={}; reports={}; end
key=sprintf('%.17g_%.17g_%.17g_%.17g_%.17g',fs,passHz,stopHz,attenuation,rippleDB);
hit=find(strcmp(keys,key),1);
if ~isempty(hit), h=coefficients{hit}; audit=reports{hit}; return; end
beta=.1102*max(0,attenuation-8.7); fc=(passHz+stopHz)/2;
order=2*ceil(max(8,(attenuation-8)/(2.285*2*pi*(stopHz-passHz)/fs))/2);
for attempt=1:12
    n=(-order/2:order/2)'; z=2*fc/fs*n; s=ones(size(z)); q=z~=0;
    s(q)=sin(pi*z(q))./(pi*z(q));
    w=besseli(0,beta*sqrt(max(0,1-(2*n/order).^2)))/besseli(0,beta);
    h=2*fc/fs*s.*w; h=h/sum(h);
    N=2^nextpow2(max(16384,16*numel(h))); H=abs(fft(h,N)); f=(0:N/2)'*fs/N; H=H(1:N/2+1);
    pass=H(f<=passHz); stop=H(f>=stopHz);
    ripple=max(abs(20*log10(max(pass,realmin)))); rejection=-20*log10(max(max(stop),realmin));
    if ripple<=rippleDB&&rejection>=attenuation, break; end
    order=2*ceil(order*1.15/2);
end
if ripple>rippleDB||rejection<attenuation
    error('tbench:AAResponse','Unable to meet the requested anti-alias response.');
end
audit=struct('inputFs',fs,'passHz',passHz,'stopHz',stopHz,'order',order, ...
    'passErrorDB',ripple,'stopAttenuationDB',rejection,'halfSupportSeconds',order/(2*fs));
if numel(keys)>=32, keys(1)=[]; coefficients(1)=[]; reports(1)=[]; end
keys{end+1}=key; coefficients{end+1}=h; reports{end+1}=audit;
end
