%%

subjNumAll = [1 3 5 10 16 17 18 20];
pFitAll = [];

for i = 1:length(subjNumAll)
    [aic, pFit] = ARCtestWvInFocusMeanZspatFilterLMSeffectPlotStack(subjNumAll(i),'LMS');
    pFitAll(i,:) = pFit;
end

%%

subjNumAll = [1 3 5 10 16 17 18 20];

pFitAll = [0.3276   -0.3714; ...
           0.0786   -0.2008; ...
           0.1051   -0.0813; ...
           0.2355   -0.4299; ...
           0.3304   -0.6937; ...
           0.2406   -0.5393; ...
           0.3305   -0.4535; ...
           0.2688   -0.3825];

optDistUnq = [1.5 2.5 3.5];
lagLead875 = [];
lagLead555 = [];

for i = 1:size(pFitAll,1)
    for j = 1:length(optDistUnq)
        lagLead875(i,j) = optDistUnq(j)*pFitAll(i,1)-pFitAll(i,2);
        lagLead555(i,j) = lagLead875(i,j)-humanWaveDefocusARC(555,875,subjNumAll(i));
    end
end

