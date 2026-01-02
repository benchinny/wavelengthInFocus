function rms = ARCfitLogGaussASFfunc(x,y,p)

% HELPER FUNCTION FOR FITTING AN ACCOMMODATION SENSITIVITY FUNCTION (ASF) 
% TO DATA FROM THE OWENS (1980) PAPER. x IS SPATIAL FREQUENCY, y IS 
% SENSITIVITY.

a1 = p(1);
m1 = p(2);
s1 = p(3);

yFit = a1.*exp(-0.5.*((log(x) - log(m1))./s1).^2);
rms = sqrt(mean((yFit-y).^2));

end