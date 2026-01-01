%% STATISTICS FOR FIGURE 3

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
foldername = fullfile(dataPath,'data','PresavedFigureData');

% PRESAVED DATA OPTIONS:
load(fullfile(foldername,'allExp1DataRGB.mat'));

optDistUnq = [1.5; 2.5; 3.5];

% REARRANGE DEFOCUS AT 555NM VECTOR
dfMean555avgColor = squeeze(mean(dfMean555all(1:10,:,:)));
% CALCULATE LAGS AND LEADS
lagLead = bsxfun(@plus,dfMean555avgColor,optDistUnq);
% COMPUTE DIFFERENCE BETWEEN LAG/LEAD AT 3.5D AND LAG/LEAD AT 1.5D
lagLeadDiff = lagLead(3,:)-lagLead(1,:);
% CALCULATE MEAN AND STD ERROR
meanLagLeadDiff = mean(lagLeadDiff);
SElagLeadDiff = 1.96.*std(lagLeadDiff)./sqrt(8);

% CALCULATING AVERAGE RESPONSE GAIN
for i = 1:size(dfMean555avgColor,2)
    % MAKE REGRESSOR VECTOR FOR INTERCEPT AND SLOPE
    X = [ones(length(optDistUnq),1) optDistUnq];
    % LINEAR REGRESS USING MATLAB BACKSLASH OPERATOR (DOCUMENTED ON
    % MATHWORKS WEBSITE)
    b(:,i) = X\(-dfMean555avgColor(:,i));
end

% MEAN AND SE OF GAIN
meanGain = mean(b(2,:));
SEgain = 1.96.*std(b(2,:))./sqrt(8);

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
