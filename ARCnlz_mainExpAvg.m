%% THIS SCRIPT CREATES AND SAVES THE 'allExp1DataRGB.mat' FILE

clear;

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

if ispc
    slash = '\';
else
    slash = '/';
end

%% LOADING DATA FROM SUBJECTS AND CONCATENATING

subjNumAll = [1 3 5 10 16 17 18 20];

wvInFocusCellAll = {}; % WAVELENGTH IN FOCUS
optDistCndAll = []; % STIMULUS DISTANCE
rgbLumNormCndAll = []; % COLOR IN RGB NORMALIZED TO MAX SINGLE-PRIMARY LUMINANCE
subjNumTag = []; % SUBJECT NUMBER TAG

for i = 1:length(subjNumAll)
    [wvInFocusCell, optDistCnd, rgbLumNormCnd] = ARCnlz_mainExpSortColor(subjNumAll(i),dataPath);
    wvInFocusCellAll = [wvInFocusCellAll wvInFocusCell];
    optDistCndAll = [optDistCndAll; -optDistCnd];
    rgbLumNormCndAll = [rgbLumNormCndAll; rgbLumNormCnd];
    subjNumTag = [subjNumTag; subjNumAll(i).*ones([size(optDistCnd,1) 1])];
end

%% SAVING

% UNCOMMENT THIS TO RECREATE THE PRESAVED EXPERIMENT DATA FILE
% save([dataPath 'data' slash 'PresavedFigureData' slash 'allExp1DataRGB.mat'],'wvInFocusCellAll','optDistCndAll','rgbLumNormCndAll','subjNumTag');
