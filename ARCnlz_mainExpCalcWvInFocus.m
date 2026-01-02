function [wvMeanAll, optDistUnq, conditionsOrderedNorm, dfMean555all] = ARCnlz_mainExpCalcWvInFocus(subjNum,dataPath)

% This function sorts data from the accommodation 
% experiment. 
%
% subjNum: subject numbers. Valid numbers: 1, 3, 5, 10, 16, 17, 18, 20
% dataPath: where the data lives
%
% wvMeanAll: wavelength-in-focus from data
% optDistUnq: unique stimulus distances
% rgbUnq: unique rgb values for different conditions
% dfMean555: estimated defocus at 555nm (to help compute lags and leads of
%            accommodation according to the usual definition)

rng(1); % FIX RANDOM SEED

% LOAD DEFOCUS VALUES, COLOR CONDITIONS, AND OPTICAL DISTANCES
[defocus875,rgbAll,optDistAll,~,~,q1,q2,q3] = ARCnlzLoadDefocusAbb(subjNum,dataPath);
% USE LCA PARAMETERS TO ESTIMATE DEFOCUS AT 555NM
defocus555 = defocus875+humanWaveDefocusParameterizedARC(555,875,q1,q2,q3);

rgbUnq = unique(rgbAll,'rows'); % UNIQUE RGB VALUES

% NORMALIZE RGB VALUES SO MAX LUMINANCE IS 1
rgbLumNorm = [];
rgbLumNorm(:,1) = (rgbUnq(:,1).^2.5)./0.2442;
rgbLumNorm(:,2) = (rgbUnq(:,2).^2.7)./0.1037;
rgbLumNorm(:,3) = (rgbUnq(:,3).^2.3)./1;
rgbLumNorm(rgbLumNorm>1) = 1;

% ORDER OF CONDITIONS FOR PLOTTING: MORE BLUE TO MORE RED WITHOUT GREEN,
% THEN MORE BLUE TO MORE RED WITH GREEN
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
                         1.00 1.00 1.00; ...
                         0.625 0   0.625];

% INDICES FOR SORTING CONDITIONS ACCORDING TO PLOT
for i = 1:size(conditionsOrderedNorm,1)
    ind(i) = find(abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
                  abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
                  abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01);
end


% MATRIX FOR STORING MEAN WAVELENGTH-IN-FOCUS VALUES
wvMeanAll = [];
dfMean555all = [];

% IMPORTANT THINGS HAPPENING IN HELPER FUNCTION: GENERATE
% PREDICTIONS OF DEFOCUS FOR EACH CONDITION
[defocus875mean, ~, optDistUnq] = ARCcalcWvInFocusHelper(defocus875,rgbAll,optDistAll);
% SORT VALUES OF DEFOCUS AT 555NM
[defocus555mean, ~, optDistUnq] = ARCcalcWvInFocusHelper(defocus555,rgbAll,optDistAll);

% SORT DATA
for i = 1:length(optDistUnq)
    % ACTUAL (NEED TO CONVERT FROM DEFOCUS TO WAVELENGTH)
    dfMean = -defocus875mean(ind,i);
    dfMean555 = -defocus555mean(ind,i);
    wvMean = humanWaveDefocusInvertParameterizedARC(875,-(dfMean+optDistUnq(i)),q1,q2,q3);
    % STORE MEAN DATA AND PREDICTIONS
    wvMeanAll(:,i) = wvMean;
    dfMean555all(:,i) = dfMean555;
end

end