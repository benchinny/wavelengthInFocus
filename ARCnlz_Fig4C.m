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
    load(fullfile(foldername,'wvMeanAndPredLMstrehl.mat'));
    aicLM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredLminusMstrehl.mat'));
    aicLminusM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredLMSstrehl.mat'));
    aicLMS = aicAll; 
elseif strcmp(modelType,'finch')
    % LOAD PRE-SAVED DATA AND MODEL FITS
    load(fullfile(foldername,'wvMeanAndPredLumdeltapass.mat'));
    aicLM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredLminusM.mat'));
    aicLminusM = aicAll;
    load(fullfile(foldername,'wvMeanAndPredLMS.mat'));
    aicLMS = aicAll;    
end

figure; 
hold on;
boxplot([0.5.*(aicLM-aicLMS)' ...
         0.5.*(aicLM-aicLminusM)']);
plot(1,0.5.*(aicLM-aicLMS)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
plot(2,0.5.*(aicLM-aicLminusM)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
ylim([-12.5 60]);
set(gca,'YTick',log([10^(-5) 10^0 10^5 10^10 10^15 10^20 10^25]));
set(gca,'YTickLabel',{'10^-5' '1' '10^5' '10^10' '10^15' '10^20' '10^25'});
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Blue-yellow vs Luminance' 'Red-green vs Luminance'});