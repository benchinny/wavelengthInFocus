%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
subjNumAll = [11 13 15 20 26 27 28 30]-10;
aicAll = [];

for i = 1:length(subjNumAll)
    aic = ARCtestWvInFocusMeanZdeltaPassPlotStack(subjNumAll(i),'LM',2,dataPath);
    aicAll(i) = aic;
end

%%

aicAll2 = [-24.2269 -33.3423 -4.7110 -14.0518 -42.0595 -44.0608 -15.9615 1.4181];

