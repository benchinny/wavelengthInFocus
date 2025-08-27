%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
subjNum = [1 3 5 10 16 17 18 20];
wvMeanAll = [];
wvPredAll = [];
aicAll = [];
dfPredPurpleAll = [];
wLMminAll = [];
wLpropMinAll = [];
pFitAll = [];

for i = 1:length(subjNum)
    [aic, pFit, wvMean, wvPred, dfPredPurple, wLMmin, wLpropMin] = ARCtestWvInFocusMeanZspatFilterLMSeffectPlotWave(subjNum(i),'LMS',dataPath);
    wvMeanAll(:,:,i) = wvMean;
    wvPredAll(:,:,i) = wvPred;
    aicAll(i) = aic;
    dfPredPurpleAll(i) = dfPredPurple;
    wLMminAll(i) = wLMmin;
    wLpropMinAll(i) = wLpropMin;
    pFitAll(:,i) = pFit;
end 

%%

optDist = [1.5 2.5 3.5];
defocusAt875fit = [];

for i = 1:size(pFitAll,2)
    for j = 1:length(optDist)
        defocusAt875fit(i,j) = -optDist(j)*pFitAll(1,i) - pFitAll(2,i);
    end
end