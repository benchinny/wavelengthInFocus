%% FIGURE 4C

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
foldername = fullfile(dataPath,'data','PresavedFigureData');

% OPTIONS FOR modelType: 'normal', 'strehl', and 'finch'
modelType = 'normal';

if strcmp(modelType,'normal')
    % LOAD PRE-SAVED DATA AND MODEL FITS
    load(fullfile(foldername,'wvMeanAndPredLM.mat'));
    aicLM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredLminusM.mat'));
    aicLminusM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredLMS.mat'));
    aicLMS = aicAll;
elseif strcmp(modelType,'strehl')
    % LOAD PRE-SAVED DATA AND MODEL FITS
    load(fullfile(foldername,'wvMeanAndPredStrehlLM.mat'));
    aicLM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredStrehlLminusM.mat'));
    aicLminusM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredStrehlLMS.mat'));
    aicLMS = aicAll; 
elseif strcmp(modelType,'finch')
    % LOAD PRE-SAVED DATA AND MODEL FITS
    load(fullfile(foldername,'wvMeanAndPredDeltaPass2.mat'));
    aicLM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredLminusM.mat'));
    aicLminusM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredLMS.mat'));
    aicLMS = aicAll;    
end

figure; 
hold on;
boxplot([(aicLM-aicLMS)' ...
          (aicLM-aicLminusM)']);
plot(1,(aicLM-aicLMS)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
plot(2,(aicLM-aicLminusM)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
ylim([-25 120]);
set(gca,'YTick',[-23.0259 0 23.0259 46.0517 69.0776 92.1034 115.1293]);
set(gca,'YTickLabel',{'10^-10' '1' '10^10' '10^20' '10^30' '10^40' '10^50'});
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Blue-yellow vs Luminance' 'Red-green vs Luminance'});