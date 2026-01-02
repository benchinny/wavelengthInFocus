function rms = ARCfitLogGaussASFfunc(x,y,p)

a1 = p(1);
m1 = p(2);
s1 = p(3);

yFit = a1.*exp(-0.5.*((log(x) - log(m1))./s1).^2);
rms = sqrt(mean((yFit-y).^2));

end