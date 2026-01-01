%% STATISTICS FOR FIGURE 3

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
foldername = fullfile(dataPath,'data','PresavedFigureData');

% PRESAVED DATA OPTIONS:
load(fullfile(foldername,'allExp1DataRGB.mat'));

optDistUnq = [1.5; 2.5; 3.5];

% REARRANGE DEFOCUS AT 555NM VECTOR (GRAB NON-CONTROL CONDITIONS)
dfMean555avgColor = squeeze(mean(dfMean555all(1:10,:,:)));
% CALCULATE LAGS AND LEADS
lagLead = bsxfun(@plus,dfMean555avgColor,optDistUnq);
% COMPUTE DIFFERENCE BETWEEN LAG/LEAD AT 3.5D AND LAG/LEAD AT 1.5D
lagLeadDiff = lagLead(3,:)-lagLead(1,:);
% CALCULATE MEAN AND 95% CIS
meanLagLeadDiff = mean(lagLeadDiff);
CI95lagLeadDiff = 1.96.*std(lagLeadDiff)./sqrt(8);

% CALCULATING AVERAGE RESPONSE GAIN
for i = 1:size(dfMean555avgColor,2)
    % MAKE REGRESSOR VECTOR FOR INTERCEPT AND SLOPE
    X = [ones(length(optDistUnq),1) optDistUnq];
    % LINEAR REGRESS USING MATLAB BACKSLASH OPERATOR (DOCUMENTED ON
    % MATHWORKS WEBSITE)
    b(:,i) = X\(-dfMean555avgColor(:,i));
end

% MEAN AND 95% CIS OF GAIN
meanGain = mean(b(2,:));
CI95gain = 1.96.*std(b(2,:))./sqrt(8);

% WAVELENGTHS IN FOCUS FOR THE CONDITIONS WITH EQUAL LUMINANCE ACROSS ALL
% PRIMARIES
wvMeanEqualLumAcross = squeeze(wvMeanAll(11,:,:));
% COMPUTE MEAN AND CI95
wvMeanMeanEqualLumAcross = mean(wvMeanEqualLumAcross,2);
wvMeanCI95equalLumAcross = 1.96.*std(wvMeanEqualLumAcross,0,2)./sqrt(8);

% GRAB CONTROL CONDITION AND ITS EQUIVALENT COMPARISON
dfMean555LumComparison = dfMean555all([3 12],:,:);
% AVERAGE ACROSS OPTICAL DISTANCE
dfMean555LumComparisonMean = squeeze(mean(dfMean555LumComparison,2));
% AVERAGE ACROSS PARTICIPANTS
dfMean555highLowMean = mean(dfMean555LumComparisonMean,2);
% 95% CIS ACROSS PARTICIPANTS
dfMean555highLowCI95 = 1.96.*std(dfMean555LumComparisonMean,0,2)./sqrt(8);
% T-TEST
[h,p,~,stats] = ttest(dfMean555LumComparisonMean(1,:)-dfMean555LumComparisonMean(2,:));

%% STATISTICS FOR FIGURE 4

clear all;

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
foldername = fullfile(dataPath,'data','PresavedFigureData');

% ALL SUBJECT NUMBERS
subjNum = [1 3 5 10 16 17 18 20];

% CELL CONTAINING PRESAVED DATA FILENAMES
presavedFilenames = {'wvMeanAndPredLM' 'wvMeanAndPredLminusM' 'wvMeanAndPredDonutx2'};

meanLagAll = []; % MEAN LAGS
CI95lagAll = []; % 95% CIS FOR LAGS
LtoMratioMeanAll = []; % MEAN L/M WEIGHT RATIO
LtoMratioCI95All = []; % 95% CIS FOR L/M WEIGHT RATIO 
StoLplusMratioMeanAll = []; % MEAN S/(L+M) WEIGHT RATIO
StoLplusMratioCI95All = []; % 95% CI S/(L+M) WEIGHT RATIO

for i = 1:length(presavedFilenames) % LOOP OVER PRESAVED FILENAMES
    load(fullfile(foldername,presavedFilenames{i}));
    % RUN HELPER FUNCTION
    if strcmp(presavedFilenames{i},'wvMeanAndPredDonutx2')
        [meanLagAll(i,:), CI95lagAll(i,:), LtoMratioMeanAll(i), LtoMratioCI95All(i), StoLplusMratioMeanAll, StoLplusMratioCI95All] = ARCnlz_Fig4statsHelper(pFitAll,wLpropMinAll,wLMminAll);
    else
        [meanLagAll(i,:), CI95lagAll(i,:), LtoMratioMeanAll(i), LtoMratioCI95All(i), StoLplusMratioMeanAll, StoLplusMratioCI95All] = ARCnlz_Fig4statsHelper(pFitAll,wLpropMinAll,[]);
    end
end
