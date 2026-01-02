function [q1,q2,q3,err] = ARC_LCAfit(wave,D)

% fits LCA function from Thibos et al. (1992) to defocus values

% SET FMINCON OPTIONS
opts             = optimset('fmincon');
opts.Algorithm   = 'active-set';
opts.LargeScale  = 'off';
% opts.UseParallel = 'never';
opts.Display     = 'off';
% opts.MaxIter     = 500;

% INITIAL CONDITIONS, UPPER AND LOWER BOUNDS
p0 = [rand*2 rand*0.99 rand*0.37];
% BOUNDS ARE NECESSARY BECAUSE MODEL IS A BIT OVER-PARAMETERIZED--3 FREE
% PARAMETERS FOR 3 DATA POINTS. WITHOUT BOUNDS, FITS WILL LOOK WEIRD (E.G.
% CURVING IN OPPOSITE DIRECTION) AND PARAMETERS VALUES WILL BE HUGELY
% VARIABLE. BUT MORE THAN HALF THE SUBJECTS WILL BUMP UP AGAINST THE BOUND
% FOR ONE VARIABLE NO MATTER HOW BOUNDS ARE SET. SO I PICKED BOUNDS THAT
% ROUGHLY BRACKETED THE VALUES IN THIBOS ET AL. (1992)
pLB = [0.01 0.01 0.01];
pUB = [2.00 0.99 0.37];

[pFit,err] = fmincon(   @(p) humanWaveDefocusObjFunc(wave,D,p),p0,[],[],[],[],pLB,pUB,[],opts);

% BEST FIT VALUES
q1 = pFit(1);
q2 = pFit(2);
q3 = pFit(3);

end