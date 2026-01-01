%% THIS SCRIPT CREATES AND SAVES THE 'allExp1DataRGB.mat' FILE

clear;

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

%% LOADING DATA FROM SUBJECTS AND CONCATENATING

subjNumAll = [1 3 5 10 16 17 18 20];
subjNumTagList = 'ABCDEFGH';

wvMeanAll = []; % WAVELENGTH IN FOCUS
optDistCndAll = []; % STIMULUS DISTANCE
conditionsOrderedNormAll = []; % COLOR IN RGB NORMALIZED TO MAX SINGLE-PRIMARY LUMINANCE
subjNumTag = []; % SUBJECT NUMBER TAG
dfMean555all = []; % ESTIMATED DEFOCUS AT 555NM (FOR LAGS AND LEADS)

for i = 1:length(subjNumAll)
    [wvMean, optDistUnq, conditionsOrderedNorm, dfMean555] = ARCnlz_mainExpCalcWvInFocus(subjNumAll(i),dataPath);
    subjNumTag = [subjNumTag; repmat(subjNumTagList(i),size(wvMean(:)))];
    wvMeanAll(:,:,i) = wvMean;
    optDistCndAll(:,i) = optDistUnq;
    conditionsOrderedNormAll(:,:,i) = conditionsOrderedNorm;
    dfMean555all(:,:,i) = dfMean555;
end

%% SAVING

% UNCOMMENT THIS TO RECREATE THE PRESAVED EXPERIMENT DATA FILE
save(fullfile(dataPath,'data','PresavedFigureData','allExp1DataRGB.mat'),'wvMeanAll','optDistCndAll','conditionsOrderedNormAll','subjNumTag','dfMean555all');
