%% STATISTICS FOR FIGURE 3

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
foldername = fullfile(dataPath,'data','PresavedFigureData');

% PRESAVED DATA OPTIONS:
load(fullfile(foldername,'allExp1DataRGB.mat'));

% REARRANGE DEFOCUS AT 555NM VECTOR
dfMean555avgColor = squeeze(mean(dfMean555all(1:10,:,:)));
% CALCULATE LAGS AND LEADS
lagLead = bsxfun(@plus,dfMean555avgColor,[1.5; 2.5; 3.5]);
% COMPUTE DIFFERENCE BETWEEN LAG/LEAD AT 3.5D AND LAG/LEAD AT 1.5D
lagLeadDiff = lagLead(3,:)-lagLead(1,:);
% CALCULATE MEAN AND STD ERROR
meanLagLeadDiff = mean(lagLeadDiff);
SElagLeadDiff = 1.96.*std(lagLeadDiff)./sqrt(8);

%% STATISTICS FOR FIGURE 4

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
foldername = fullfile(dataPath,'data','PresavedFigureData');

% ALL SUBJECT NUMBERS
subjNum = [1 3 5 10 16 17 18 20];

% PRESAVED DATA OPTIONS:
% wvMeanAndPredLminusM: RED-GREEN PREDICTIONS
% wvMeanAndPredDonutx2: BLUE-YELLOW PREDICTIONS
% wvMeanAndPredLM: LUMINANCE PREDICTIONS
load(fullfile(foldername,'wvMeanAndPredLM.mat'));

symbDist = 'sod'; % SYMBOLS FOR PLOTTING
% ORDER CONDITIONS FOR PLOTTING
conditionsOrderedNorm = [0.25 0.00 1.00; ...
                         0.50 0.00 1.00; ...
                         1.00 0.00 1.00; ...
                         1.00 0.00 0.50; ...
                         1.00 0.00 0.25; ...
                         0.25 0.50 1.00; ...
                         0.50 0.50 1.00; ...
                         1.00 0.50 1.00; ...
                         1.00 0.50 0.50; ...
                         1.00 0.50 0.25; ...
                         1.00 1.00 1.00];
