%% STATISTICS FOR FIGURE 3

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
foldername = fullfile(dataPath,'data','PresavedFigureData');

% PRESAVED DATA OPTIONS:
load(fullfile(foldername,'allExp1DataRGB.mat'));

optDistUnq = [1.5; 2.5; 3.5]; % UNIQUE OPTICAL DISTANCES

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

% DISPLAY STATEMENTS BELOW

display('----------------');
display(['mean lag difference between 3.5D and 1.5D = ' num2str(meanLagLeadDiff,3) ' ± ' num2str(CI95lagLeadDiff,3) newline]);

display(['mean gain = ' num2str(meanGain,3) ' ± ' num2str(CI95gain,3) newline]);

display(['wavelength in focus when luminance equal across all display primaries = ']);
display([num2str(wvMeanMeanEqualLumAcross(1),3) 'nm ± ' num2str(wvMeanCI95equalLumAcross(1),3) 'nm ' 'when stim distance = ' num2str(optDistUnq(1)) 'D']);
display([num2str(wvMeanMeanEqualLumAcross(2),3) 'nm ± ' num2str(wvMeanCI95equalLumAcross(2),3) 'nm ' 'when stim distance = ' num2str(optDistUnq(2)) 'D']);
display([num2str(wvMeanMeanEqualLumAcross(3),3) 'nm ± ' num2str(wvMeanCI95equalLumAcross(3),3) 'nm ' 'when stim distance = ' num2str(optDistUnq(3)) 'D' newline]);

display(['mean accommodation for lower luminance = ' num2str(-dfMean555highLowMean(2),3) ' ± ' num2str(dfMean555highLowCI95(2),3)]);
display(['mean accommodation for higher luminance = ' num2str(-dfMean555highLowMean(1),3) ' ± ' num2str(dfMean555highLowCI95(1),3) newline]);
display(['t(' num2str(stats.df) ')= ' num2str(stats.tstat,3) ', p = ' num2str(p,3) newline]);

%% STATISTICS FOR FIGURE 4

% ALL SUBJECT NUMBERS
subjNum = [1 3 5 10 16 17 18 20];
optDistUnq = [1.5; 2.5; 3.5]; % UNIQUE OPTICAL DISTANCES

% CELL CONTAINING PRESAVED DATA FILENAMES
presavedFilenames = {'wvMeanAndPredLM' 'wvMeanAndPredLminusM' 'wvMeanAndPredLMS'};

meanLagAll = []; % MEAN LAGS
CI95lagAll = []; % 95% CIS FOR LAGS
LtoMratioMeanAll = []; % MEAN L/M WEIGHT RATIO
LtoMratioCI95All = []; % 95% CIS FOR L/M WEIGHT RATIO 
StoLplusMratioMeanAll = []; % MEAN S/(L+M) WEIGHT RATIO
StoLplusMratioCI95All = []; % 95% CI S/(L+M) WEIGHT RATIO

for i = 1:length(presavedFilenames) % LOOP OVER PRESAVED FILENAMES
    load(fullfile(foldername,presavedFilenames{i}));
    % RUN HELPER FUNCTION
    if strcmp(presavedFilenames{i},'wvMeanAndPredLMS')
        [meanLagAll(i,:), CI95lagAll(i,:), LtoMratioMeanAll(i), LtoMratioCI95All(i), StoLplusMratioMeanAll, StoLplusMratioCI95All] = ARCnlz_Fig4statsHelper(pFitAll,wLpropMinAll,wLMminAll,wSall);
    else
        [meanLagAll(i,:), CI95lagAll(i,:), LtoMratioMeanAll(i), LtoMratioCI95All(i), StoLplusMratioMeanAll, StoLplusMratioCI95All] = ARCnlz_Fig4statsHelper(pFitAll,wLpropMinAll,[],[]);
    end
end

% DISPLAY STATEMENTS BELOW

display('----------------');
display(['Mean lag parameters: ']);
display(['L+M model: ' num2str(meanLagAll(1,1),3) ' ± ' num2str(CI95lagAll(1,1),3) ' for ' num2str(optDistUnq(1)) 'D']);
display(['L+M model: ' num2str(meanLagAll(1,2),3) ' ± ' num2str(CI95lagAll(1,2),3) ' for ' num2str(optDistUnq(2)) 'D']);
display(['L+M model: ' num2str(meanLagAll(1,3),3) ' ± ' num2str(CI95lagAll(1,3),3) ' for ' num2str(optDistUnq(3)) 'D']);
display(['L-M model: ' num2str(meanLagAll(2,1),3) ' ± ' num2str(CI95lagAll(2,1),3) ' for ' num2str(optDistUnq(1)) 'D']);
display(['L-M model: ' num2str(meanLagAll(2,2),3) ' ± ' num2str(CI95lagAll(2,2),3) ' for ' num2str(optDistUnq(2)) 'D']);
display(['L-M model: ' num2str(meanLagAll(2,3),3) ' ± ' num2str(CI95lagAll(2,3),3) ' for ' num2str(optDistUnq(3)) 'D']);
display(['(L+M)-S model: ' num2str(meanLagAll(3,1),3) ' ± ' num2str(CI95lagAll(3,1),3) ' for ' num2str(optDistUnq(1)) 'D']);
display(['(L+M)-S model: ' num2str(meanLagAll(3,2),3) ' ± ' num2str(CI95lagAll(3,2),3) ' for ' num2str(optDistUnq(2)) 'D']);
display(['(L+M)-S model: ' num2str(meanLagAll(3,3),3) ' ± ' num2str(CI95lagAll(3,3),3) ' for ' num2str(optDistUnq(3)) 'D' newline]);
display(['L/M weight ratio: ']);
display(['L+M model: ' num2str(LtoMratioMeanAll(1),3) ' ± ' num2str(LtoMratioCI95All(1),3)]);
display(['L-M model: ' num2str(LtoMratioMeanAll(2),3) ' ± ' num2str(LtoMratioCI95All(2),3)]);
display(['(L+M)-S model: ' num2str(LtoMratioMeanAll(3),3) ' ± ' num2str(LtoMratioCI95All(3),3) newline]);

display(['S/(L+M) weight ratio: ']);
display(['(L+M)-S model: ' num2str(StoLplusMratioMeanAll,3) ' ± ' num2str(StoLplusMratioCI95All,3) newline]);