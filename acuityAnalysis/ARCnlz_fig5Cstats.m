%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
subjNumAll = [1 3 5 10 16 17 18 20];

pNull = [];
for i = 1:length(subjNumAll)
    [unqFocDst,PC,PCci,dprime,dprimeCI,PCfit,dprimeFitAll, ...
    PCfitSupport,bestDist,bestDistCI,PCboots] = ...
    ARCacuAnalysisSubjective(subjNumAll(i),0,dataPath);   
    bPeakNot0 = [];
    for j = 1:size(PCboots,2)
        PCtmp = PCboots(:,j);
        [PCmaxTmp,indPeak] = max(PCtmp);
        bPeakNot0(j) = PCmaxTmp>PCtmp(5);
        epsilonPC = 0.99;
        PCforDP = PCtmp;
        PCforDP(PCforDP>epsilonPC) = epsilonPC;
        dprime = 2*norminv(PCforDP);
    end
    pNull(i) = 1-mean(bPeakNot0);
end