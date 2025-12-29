function [stdEst,LL] = ARCfitStdGauss(res)

% SET FMINCON OPTIONS
opts             = optimset('fmincon');
opts.Algorithm   = 'active-set';
opts.LargeScale  = 'off';
% opts.UseParallel = 'never';
opts.Display     = 'off';
% opts.MaxIter     = 500;

p0 = rand;
[stdEst,LL] = fmincon(   @(p) ARCfitStdGaussFunc(res,p),p0,[],[],[],[],[],[],[],opts);

end