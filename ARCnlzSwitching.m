function [d, wS, rbThreshold, mse] = ARCnlzSwitching(deltaA,deltaRB1,deltaRB2,deltaS,dMin)

% SET FMINCON OPTIONS
opts             = optimset('fmincon');
opts.Algorithm   = 'active-set';
opts.LargeScale  = 'off';
% opts.UseParallel = 'never';
opts.Display     = 'iter';
opts.MaxIter     = 500;

ub = [1.5 2    0.5*1];
lb = [dMin  -0.5 -0.5*1];
p0 = [rand*1.5 rand*2 (rand-0.5)*1];
[pFit,mse] = fmincon(   @(p) ARCnlzSwitchingFunc(deltaA,deltaRB1,deltaRB2,deltaS,p),p0,[],[],[],[],lb,ub,[],opts);

d = pFit(1);
wS = pFit(2);
rbThreshold = pFit(3);

end