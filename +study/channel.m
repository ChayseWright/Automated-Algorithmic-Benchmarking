function c = channel(ctx,m)
%CHANNEL Single-muscle view with independent reference and reduced raw branch.
c=tbench.channel(ctx,m);
c.reducedRaw=ctx.reducedRaw(:,m); c.reducedRawTime=ctx.reducedRawTime;
c.reducedRawFs=ctx.reducedRawFs; c.observerGuardSeconds=ctx.observerGuardSeconds;
end
