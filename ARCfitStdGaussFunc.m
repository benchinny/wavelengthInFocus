function LL = ARCfitStdGaussFunc(res,p)

LL = -sum(log(normpdf(res,0,p)));

end