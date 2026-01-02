function [pFit,err] = ARCfitLagLead(x,y,d,bFITLINE,objFunc)

% function for fitting lags and lead parameters in the context of the
% signal quality model shown in Figure 4. What's going on conceptually is
% described in the methods section of the manuscript.
%
% x: predicted accommodation from model (in diopters)
% y: actual accommodation
% d: accommodative demand
% bFITLINE: whether or not to use assumption that lags and leads grow
%           linearly
% objFunc : whether or not to use RMSE or negative log-likelihood as the
%           objective function (basically give same result)
%
% pFit: lag / lead parameters
% err : error function

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
    [pFit,err] = fmincon(   @(p) ARCfitLagLeadFuncLin(x,y,d,p,objFunc),p0,[],[],[],[],lb,ub,[],opts);
else
    p0 = rand([1 length(unique(d))]);
    lb = -20.*ones([1 length(unique(d))]);
    [pFit,err] = fmincon(   @(p) ARCfitLagLeadFunc(x,y,d,p,objFunc),p0,[],[],[],[],lb,[],[],opts);    
end

end