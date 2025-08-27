%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
subjNumAll = [11 13 15 20 26 27 28 30];
aicAll = [];

for i = 1:length(subjNumAll)
    aic = ARCtestWvInFocusMeanZdeltaPassPlotStack(subjNumAll(i),'LM',2,dataPath);
    aicAll(i) = aic;
end