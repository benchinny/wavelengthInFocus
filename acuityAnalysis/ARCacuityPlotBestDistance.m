%% MAKE FIGURE 5D

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
% OPTIONS:
% acuityModelingPrediction: FOR BEST-FITTING COLOR-OPPONENT MODEL
% acuityModelingPredictionLum: FOR BEST-FITTING LUMINANCE MODEL
fileStr = 'acuityModelingPrediction';

subjNumAll = [1 3 5 10 16 17 18 20];
peakLocModelPredictionAll = []; % PREDICTED 'BEST ACUITY' DISTANCES FROM MODEL
subjNumInclude = [1:4 5:8]; % FOR IF YOU WANT TO LEAVE OUT SUBJECTS FROM ANALYSIS
plotColorCorr = [1 1 1]; % PLOTTING COLOR
plotSymb = '>^<vsdho'; % SYMBOLS FOR PLOTTING INDIVIDUALS

if ispc
    slash = '\';
else
    slash = '/';
end
dataFolder = [dataPath 'data' slash 'acuityModeling' slash];

shiftValBestFitAll = []; % IF USING DEPTH-OF-FOCUS PARAMETER
for j = 1:length(subjNumAll)

    subjNum = subjNumAll(j);
    
    % LOAD PRE-GENERATED D-PRIME DATA FOR ACUITY TASK
    load([dataFolder fileStr 'S' num2str(subjNum) '.mat']);

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
    
% HARD-CODED ACTUAL PEAK LOCATIONS AND UPPER/LOWER BOUNDS
peakLocActualAll = [1.8635 1.9435 2.6335 2.1235 2.6435 1.9035 1.8935 2.3235];
peakLocLBall = [1.6835 1.8135 1.5935 1.9435 2.2535 1.7935 1.5035 2.2635];
peakLocUBall = [3.0935 3.1835 3.5635 3.1035 2.7735 2.0535 2.4535 2.4235];

figure;
hold on;
for i = 1:length(peakLocModelPredictionAll)
   plot(peakLocModelPredictionAll(i),peakLocActualAll(subjNumInclude(i)),['k' plotSymb(i)],'MarkerFaceColor',plotColorCorr,'MarkerSize',15);
end
RMSE = sqrt(mean((peakLocModelPredictionAll-peakLocActualAll(subjNumInclude)).^2));
axis square;
set(gca,'FontSize',15);
set(gca,'Box','on');
set(gca,'XTick',1.6:0.2:2.8);
set(gca,'YTick',1.6:0.2:2.8);
xlim([1.5 2.8]);
ylim([1.5 2.8]);
plot([1.5 2.8],[1.5 2.8],'k--','LineWidth',1);
xlabel(['Predicted Peak Location (D)']);
ylabel(['Actual Peak Location (D)']);
title(['Correlation = ' num2str(corr(peakLocModelPredictionAll',peakLocActualAll(subjNumInclude)'),3)]);

%% T-TEST ERROR VALUES

errorLum = [-0.2035 -0.0335 -0.4135 0.0965 ...
            -0.9935 0.1165 0.1665 -0.5135];
errorCO = [-0.3035 -0.1335 -0.5435 -0.1735 ...
           -0.7135 0.0765 -0.0135 -0.4835];

[h,p] = ttest(errorCO-errorLum);