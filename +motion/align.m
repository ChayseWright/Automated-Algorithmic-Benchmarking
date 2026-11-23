function [y,good] = align(time,x,target,m)
%ALIGN Reduce before interpolation when motion has a lower sample rate.
fs=1/median(diff(time)); query=target.time; fsTarget=target.fs;
if fsTarget<fs*(1-1e-7)
    pass=max(1.1*m.searchBand(2),min(30,.36*fsTarget)); stop=.45*fsTarget;
    [x,t]=tbench.preprocess.rateChange(x,fs,fsTarget,pass,stop,m.reduction);
    time=time(1)+t;
end
[y,good]=tbench.motion.sample(time,x,query,m.reduction.sincHalfLength);
end
