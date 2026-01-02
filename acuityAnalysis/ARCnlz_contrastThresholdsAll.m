function thresholdsAll = ARCnlz_contrastThresholdsAll(dataPath,bPLOT)

% this is a wrapper function for ARCnlz_contrastThresholds, which
% calculates thresholds for each subject in the calibration task for the
% LCA measurement and acuity tasks.
%
% dataPath: path to data. e.g. 
%           dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\'
% bPLOT   : plot or not
%
% thresholdsAll: matrix containing contrast thresholds for each of the 8 
%                subjects for the red, green, blue, and purple gratings
%                respectively 

% ALL SUBJECTS
subjNumAll = [1 3 5 10 16 17 18 20];
% INITIALIZE MATRIX FOR STORING THRESHOLDS
thresholdsAll = [];

for i = 1:length(subjNumAll) % LOOP OVER EACH SUBJECT
    % CALCULATE AND STORE THRESHOLDS
    thresholds = ARCnlz_contrastThresholds(subjNumAll(i),bPLOT,dataPath);
    thresholdsAll(i,:) = thresholds;
end

save(fullfile(dataPath,'data','PresavedFigureData','thresholdsAll.mat'),'thresholdsAll');

end