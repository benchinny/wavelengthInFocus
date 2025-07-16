function val = ARCfitLagLeadFunc(x,y,d,p,objFunc)

dUnq = unique(d);
if length(dUnq)~=length(p)
    error('ARCfitLagLeadFunc: fewer lag/lead parameters than unique stimulus distances!');
end

for i = 1:length(dUnq)
    ind = abs(d-dUnq(i))<0.0001;
    y(ind) = y(ind)+p(i);
end

if strcmp(objFunc,'RMS')
    val = sqrt(mean((y-x).^2));
elseif strcmp(objFunc,'NLL')
    estStd = std(y-x);
    val = -sum(log(normpdf(y-x,0,estStd)));
else
    error('ARCfitLagLeadFuncLin: invalid objective function! Either RMS or NLL');
end

end