function [pFit,rms] = ARCfitLogGaussASF(x,y)

% FITS AN ACCOMMODATION SENSITIVITY FUNCTION (ASF) TO DATA FROM THE OWENS
% (1980) PAPER. x IS SPATIAL FREQUENCY, y IS SENSITIVITY.

% SET FMINCON OPTIONS
opts             = optimset('fmincon');
opts.Algorithm   = 'active-set';
opts.LargeScale  = 'off';
% opts.UseParallel = 'never';
opts.Display     = 'off';
opts.MaxIter     = 500;

p0 = rand([1 3]);
lb = [0.01 0.01 0.01];
[pFit,rms] = fmincon(   @(p) ARCfitLogGaussASFfunc(x,y,p),p0,[],[],[],[],lb,[],[],opts);

end