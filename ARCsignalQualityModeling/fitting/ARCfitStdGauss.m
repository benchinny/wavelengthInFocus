function [stdEst,LL] = ARCfitStdGauss(res)

% function for estimating standard deviation of data in a manner consistent
% with the way MATLAB does it to calculate AIC internally. Basically, finds
% the standard deviation that maximizes the log-likelihood of a Gaussian
% fit to some data, assuming mean 0.

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