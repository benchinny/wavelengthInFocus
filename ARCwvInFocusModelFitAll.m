function ARCwvInFocusModelFitAll(subjNumAll,modelType,sigQualType,dataPath)

% wrapper function for fitting main model to data
%
% subjNumAll: all subject numbers. Valid numbers: 1, 3, 5, 10, 16, 17, 18,
%                                                 20
% modelType: which color model to use
%           'LMS': blue-yellow model.
%           'LM': luminance model
%           'LminusM': red-green model
% sigQualType: which signal quality metric to use
%           'xcorr': cross-correlation metric.
%           'strehl': Strehl ratio
% dataPath: path to folder with all data. Example below.
%           dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

for i = 1:length(subjNumAll)
    ARCwvInFocusModelFit(subjNumAll(i),modelType,sigQualType,dataPath);
end

end


