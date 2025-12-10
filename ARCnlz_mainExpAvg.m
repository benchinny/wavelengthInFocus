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

wvInFocusCellAll = {};
defocusAt550cellAll = {};
defocusAt875cellAll = {};
optDistCndAll = [];
rgbLumNormCndAll = [];
subjNumTag = [];

for i = 1:length(subjNumAll)
    [wvInFocusCell, defocusAt550cell, defocusAt875cell, optDistCnd, rgbLumNormCnd] = ARCnlz_mainExpSortColor(subjNumAll(i)+10,dataPath);
    wvInFocusCellAll = [wvInFocusCellAll wvInFocusCell];
    defocusAt550cellAll = [defocusAt550cellAll defocusAt550cell];
    defocusAt875cellAll = [defocusAt875cellAll defocusAt875cell];
    optDistCndAll = [optDistCndAll; -optDistCnd];
    rgbLumNormCndAll = [rgbLumNormCndAll; rgbLumNormCnd];
    subjNumTag = [subjNumTag; subjNumAll(i).*ones([size(optDistCnd,1) 1])];
end

%% SAVING

% UNCOMMENT THIS TO RECREATE THE PRESAVED EXPERIMENT DATA FILE
% save([dataPath 'data' slash 'PresavedFigureData' slash 'allExp1DataRGB.mat'],'wvInFocusCellAll','defocusAt550cellAll','defocusAt875cellAll','optDistCndAll','rgbLumNormCndAll','subjNumTag');
