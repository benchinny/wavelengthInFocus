function [RMSE, corrPred, p, peakLocModelPredictionAll, peakLocActualAll] = ARCnlz_Fig5Dhelper(dataPath,fileStr)

% helper function for making Figure 5D and computing stats
%
% dataPath: folder where the 'data' folder lives
% fileStr: string of pre-saved data in the 'acuityModeling' folder

% OPTIONS:
% acuityModelingPrediction: FOR BEST-FITTING COLOR-OPPONENT MODEL
% acuityModelingPredictionLM: FOR BEST-FITTING LUMINANCE MODEL

subjNumAll = [1 3 5 10 16 17 18 20]; % LIST OF ALL SUBJECT NUMBERS
peakLocModelPredictionAll = []; % PREDICTED 'BEST ACUITY' DISTANCES FROM MODEL
subjNumInclude = [1:4 5:8]; % FOR IF YOU WANT TO LEAVE OUT SUBJECTS FROM ANALYSIS

dataFolder = fullfile(dataPath,'data','acuityModeling');

shiftValBestFitAll = []; % IF USING DEPTH-OF-FOCUS PARAMETER
for j = 1:length(subjNumAll)

    subjNum = subjNumAll(j);
    
    % LOAD PRE-GENERATED D-PRIME DATA FOR ACUITY TASK
    load(fullfile(dataFolder,[fileStr 'S' num2str(subjNum) '.mat']));

    scaleFac = 0.816; % CONVERT OPTOTUNE VALUE IN BVAMS TO DEFOCUS AT EYE
    
    % 'DEPTH-OF-FOCUS' PARAMETER (SET TO 0 CURRENTLY)
    shiftVals = 0;
    
    % SHIFT PREDICTED D-PRIME CURVES AND REGRESS AGAINST ACTUAL DATA TO
    % FIT 'DEPTH-OF-FOCUS' PARAMETER
    for i = 1:length(shiftVals)
        % INTERPOLATE MODEL PREDICTIONS TO BE AT SAME X-VALUES AS
        % PARTICIPANTS' DATA        
        dprime2regressTmp = interp1(defocusForStim+modelPrediction875nmPurpleAt2pt5+shiftVals(i),dprimeMetric,2.5+unqFocDst.*scaleFac);
        % CALCULATE BEST-FITTING SCALE FACTOR
        dprimeScaleTmp(i) = dprime2regressTmp\dprime';
        % STORE ERROR
        errorDP(i) = sqrt(mean((dprimeScaleTmp(i).*dprime2regressTmp-dprime').^2));
    end
    % FIND SHIFT YIELDING BEST FOCUS
    [~,indMinShift] = min(errorDP);
    shiftValBestFit = shiftVals(indMinShift);
    shiftValBestFitAll(j) = shiftValBestFit;
    
    % DETERMINE PREDICTED PEAK PERFORMANCE AFTER FITTING DEPTH-OF-FOCUS
    % PARAMETER    
    stimDistanceSmp = 1.2:0.01:3.8;
    dprimeMetricSmooth = interp1(defocusForStim+modelPrediction875nmPurpleAt2pt5+shiftVals(indMinShift),dprimeMetric,stimDistanceSmp,'spline');
    [~,indPeak] = max(dprimeMetricSmooth);
    peakLocModelPrediction = stimDistanceSmp(indPeak);
    peakLocModelPredictionAll(j) = peakLocModelPrediction;
end
   
% GET EMPIRICALLY-ESTIMATED BEST DISTANCES
peakLocActualAll = [];
for i = 1:length(subjNumAll)
    [~,~,~,~,~,~,~,~,bestDist,bestDistCI,~] = ARCacuityAnalyzeDataOnly(subjNumAll(i),0,dataPath);
    peakLocActualAll(i) = 2.5+bestDist;
end

% CALCULATE RMSE
RMSE = sqrt(mean((peakLocModelPredictionAll-peakLocActualAll(subjNumInclude)).^2));
% CALCULATE CORRELATION
[corrPred,p] = corr(peakLocModelPredictionAll',peakLocActualAll(subjNumInclude)');
