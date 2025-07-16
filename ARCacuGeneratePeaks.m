%%

subjNumAll = [1 3 5 10 16 17 18 20];
bestDistAll = [];
bestDistCIall = [];

for i = 1:length(subjNumAll)
    [unqFocDst,PC,PCci,dprime,dprimeCI,PCfit, ...
    dprimeFitAll,PCfitSupport,bestDist,bestDistCI] = ...
    ARCacuAnalysisSubjective(subjNumAll(i),0);
    bestDistAll(i,:) = bestDist;
    bestDistCIall(i,:) = bestDistCI;
end