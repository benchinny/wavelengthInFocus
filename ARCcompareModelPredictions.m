%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];

load([foldername 'wvMeanAndPredDonutx2.mat']);
dfPredPurpleAllBY = dfPredPurpleAll;
aicBY = aicAll;
load([foldername 'wvMeanAndPredLminusM.mat']);
dfPredPurpleAllRG = dfPredPurpleAll;
aicRG = aicAll;
load([foldername 'wvMeanAndPredLM.mat']);
dfPredPurpleAllLum = dfPredPurpleAll;
aicLum = aicAll;

figure; 
set(gcf,'Position',[1 1 1132 539]);
subplot(1,3,1);
plot(-dfPredPurpleAllBY,-dfPredPurpleAllRG,'ko','MarkerSize',10, ...
    'MarkerFaceColor','w'); 
axis square; 
formatFigure('Blue-Yellow','Red-Green',['\rho = ' num2str(corr(dfPredPurpleAllBY',dfPredPurpleAllRG'),2)]);
xlim([1 2]);
ylim([1 2]);
subplot(1,3,2);
plot(-dfPredPurpleAllLum,-dfPredPurpleAllBY,'ko','MarkerSize',10, ...
    'MarkerFaceColor','w'); 
axis square; 
formatFigure('Luminance','Blue-yellow',['\rho = ' num2str(corr(dfPredPurpleAllLum',dfPredPurpleAllBY'),2)]);
xlim([1 2]);
ylim([1 2]);
subplot(1,3,3);
plot(-dfPredPurpleAllLum,-dfPredPurpleAllRG,'ko','MarkerSize',10, ...
    'MarkerFaceColor','w'); 
axis square; 
formatFigure('Luminance','Red-green',['\rho = ' num2str(corr(dfPredPurpleAllLum',dfPredPurpleAllRG'),2)]);
xlim([1 2]);
ylim([1 2]);

