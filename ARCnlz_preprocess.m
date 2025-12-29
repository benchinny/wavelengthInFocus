%% SCRIPT FOR PREPROCESSING--GRATING CONTRAST THRESHOLDS AND LCA

% This script performs preprocessing steps and saves out data necessary for
% running the rest of the analyses shown in the Chin et al. (2026) paper. 

% MODIFY THIS PATH ON YOUR COMPUTER--IT SHOULD POINT TO WHEREVER THE 'data'
% FOLDER YOU DOWNLOADED FROM ZENODO IS
dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
% WHETHER OR NOT TO PLOT--CHANGE TO true IF YOU WANT TO SEE THE PLOTS
bPLOT = true;

% GENERATE EACH SUBJECT'S CONTRAST THRESHOLDS FOR DIFFERENT COLORED GRATINGS
% (ONLY RELEVANT TO ACUITY MODEL)
thresholdsAll = ARCnlz_contrastThresholdsAll(dataPath,bPLOT);

% FIT LCA PARAMETERS FOR EACH SUBJECT (AND MAKES FIGURE 7E FROM PAPER)
[q1bestAll, q2bestAll, q3bestAll] = ARCnlzLCAall(dataPath,bPLOT);