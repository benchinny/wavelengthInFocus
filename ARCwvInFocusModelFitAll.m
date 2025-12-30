%% WRAPPER SCRIPT FOR FITTING MAIN MODEL TO DATA

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
% OPTIONS: 'LMS': blue-yellow model.
%           'LM': luminance model
%          'LminusM': red-green model
modelType = 'LMS';
% OPTIONS: 'xcorr': cross-correlation metric.
%           'strehl': Strehl ratio
sigQualType = 'xcorr';

ARCwvInFocusModelFit(1,modelType,sigQualType,dataPath);

ARCwvInFocusModelFit(3,modelType,sigQualType,dataPath);

ARCwvInFocusModelFit(5,modelType,sigQualType,dataPath);

ARCwvInFocusModelFit(10,modelType,sigQualType,dataPath);

ARCwvInFocusModelFit(16,modelType,sigQualType,dataPath);

ARCwvInFocusModelFit(17,modelType,sigQualType,dataPath);

ARCwvInFocusModelFit(18,modelType,sigQualType,dataPath);

ARCwvInFocusModelFit(20,modelType,sigQualType,dataPath);