%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
subjNum = [1 3 5 10 16 17 18 20];
wvMeanAll = [];
wvPredAll = [];
aicAll = [];
dfPredPurpleAll = [];
wLMminAll = [];
wLpropMinAll = [];

for i = 1:length(subjNum)
    [aic, pFit, wvMean, wvPred, dfPredPurple, wLMmin, wLpropMin] = ARCtestWvInFocusMeanZspatFilterLMSeffectPlotWave(subjNum(i),'LMS',dataPath);
    wvMeanAll(:,:,i) = wvMean;
    wvPredAll(:,:,i) = wvPred;
    aicAll(i) = aic;
    dfPredPurpleAll(i) = dfPredPurple;
    wLMminAll(i) = wLMmin;
    wLpropMinAll(i) = wLpropMin;
end 
