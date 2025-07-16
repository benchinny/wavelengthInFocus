function [pFit,rms] = ARCfitLagLead(x,y,d,bFITLINE,objFunc)

% SET FMINCON OPTIONS
opts             = optimset('fmincon');
opts.Algorithm   = 'active-set';
opts.LargeScale  = 'off';
% opts.UseParallel = 'never';
opts.Display     = 'off';
opts.MaxIter     = 500;

if bFITLINE
    p0 = rand([1 2]);
    lb = [-5 -5];
    ub = [5 5];
    [pFit,rms] = fmincon(   @(p) ARCfitLagLeadFuncLin(x,y,d,p,objFunc),p0,[],[],[],[],lb,ub,[],opts);
else
    p0 = rand([1 length(unique(d))]);
    lb = -20.*ones([1 length(unique(d))]);
    [pFit,rms] = fmincon(   @(p) ARCfitLagLeadFunc(x,y,d,p,objFunc),p0,[],[],[],[],lb,[],[],opts);    
end

end