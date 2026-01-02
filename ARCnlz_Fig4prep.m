%% GENERATE MODEL PREDICTIONS TOGETHER WITH ACTUAL DATA

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
modelType = 'LMS';
sigQualType = 'xcorr';
bSAVE = true;

savePath = fullfile(dataPath,'data','PresavedFigureData/');
subjNum = [1 3 5 10 16 17 18 20];
wvMeanAll = []; % MEASURED WAVELENGTH IN FOCUS
wvPredAll = []; % PREDICTED WAVELENGTH IN FOCUS
aicAll = []; % AIC VALUES FOR MODEL EVALUATION
dfPredPurpleAll = []; % PREDICTION FOR THE ACUITY TASK (USED FOR FIG 5)
wLMminAll = []; % FIT RATIO OF (L+M) TO S CONE WEIGHTS
wLpropMinAll = []; % RATIO OF L TO M CONE WEIGHTS
pFitAll = []; % LAG AND LEAD PARAMETERS
wSall = []; % FOR STORING WEIGHTS ON S-CONE IMAGE

if strcmp(modelType,'LMS')
    savename = 'wvMeanAndPredLMS';
elseif strcmp(modelType,'LminusM')
    savename = 'wvMeanAndPredLminusM';
elseif strcmp(modelType,'LM')
    savename = 'wvMeanAndPredLM';
elseif strcmp(modelType,'Lum')
    savename = 'wvMeanAndPredLum';
end

if strcmp(sigQualType,'strehl')
    savename = [savename 'strehl'];
elseif strcmp(sigQualType,'deltapass')
    savename = [savename 'deltapass'];    
end

for i = 1:length(subjNum)
    [aic, pFit, wvMean, wvPred, dfPredPurple, wLMmin, wLpropMin, wS] = ARCwvInFocusModelSort(subjNum(i),modelType,sigQualType,dataPath);
    wvMeanAll(:,:,i) = wvMean;
    wvPredAll(:,:,i) = wvPred;
    aicAll(i) = aic;
    dfPredPurpleAll(i) = dfPredPurple;
    wLMminAll(i) = wLMmin;
    wLpropMinAll(i) = wLpropMin;
    pFitAll(:,i) = pFit;
    wSall(i) = wS;
end 

if bSAVE
    save(fullfile(savePath,savename),'wvMeanAll', ...
        'wvPredAll','aicAll','dfPredPurpleAll', ...
        'wLMminAll','wLpropMinAll','pFitAll');
end