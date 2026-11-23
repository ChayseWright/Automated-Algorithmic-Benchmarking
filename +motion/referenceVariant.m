function t = referenceVariant(t,name,ref,m)
%REFERENCEVARIANT Prespecified positive-state/frequency sensitivity witnesses.
% All variants retain the common signed-coordinate power target. Standalone
% oscillator negatives have no amplitude-limit certificate and stay unknown.
switch name
    case 'consensus', return;
    case 'multitaper_only', t.state=t.multitaperState;
    case 'oscillator_only'
        t.state=t.oscillatorState; t.frequencyTrack=t.secondary.frequency;
        t.negativeReferenceCalibrated=false;
    otherwise, error('tbench:MotionReferenceVariant','Unknown reference sensitivity witness.');
end
t.frequency=NaN; if any(t.state==1), t.frequency=median(t.frequencyTrack(t.state==1)); end
t.positiveSeconds=sum(t.state==1)/t.fs; t.negativeSeconds=sum(t.state==-1)/t.fs;
t.indeterminateSeconds=sum(t.valid&t.state==0)/t.fs;
t.stateCoverage=sum(t.state~=0)/max(sum(t.valid),1);
p=t.powerSigned; p(isfinite(p)&p<0)=0;
t.events=tbench.episodes(t.state==1,t.time,t.frequencyTrack,p,t.valid,0);
for k=1:numel(t.events), t.events(k).boundaryUncertaintySeconds=m.windowSeconds/2+ref.clock.uncertaintySeconds; end
end
