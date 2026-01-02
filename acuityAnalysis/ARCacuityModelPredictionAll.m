%% SCRIPT FOR PRE-GENERATING DATA FOR FIGURE 5

% This script pre-generates mode predictions for Figure 5. On a computer 
% with 8 cores, the script will take around 4-5 hours to run. 

% WHERE ALL THE DATA IS LOCATED (CHANGE IF NECESSARY)
dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
% LIST OF SUBJECT NUMBERS
subjNumAll = [1 3 5 10 16 17 18 20];

% GENERATE MODEL PREDICTIONS FOR BEST CHROMATIC MODEL
for i = 1:length(subjNumAll)
    ARCacuityModelPrediction(subjNumAll(i),'Chrom',dataPath);
end

% GENERATE MODEL PREDICTIONS FOR L+M MODEL
for i = 1:length(subjNumAll)
    ARCacuityModelPrediction(subjNumAll(i),'LM',dataPath);
end