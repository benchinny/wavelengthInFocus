function LL = ARCfitStdGaussFunc(res,p)

% helper function for estimating standard deviation of data in a manner 
% consistent with the way MATLAB does it to calculate AIC internally. 
% Basically, finds the standard deviation that maximizes the 
% log-likelihood of a Gaussian fit to some data, assuming mean 0.

LL = -sum(log(normpdf(res,0,p)));

end