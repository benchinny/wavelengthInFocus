%% RECREATES FIGURE 5B-STYLE PLOTS FOR ALL SUBJECTS

% LIST OF ALL SUBJECT NUMBERS
subjNumAll = [1 3 5 10 16 17 18 20];
peakLocModelPredictionAll = []; % MODEL PREDICTIONS OF BEST ACUITY DISTANCE
dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
% OPTIONS:
% acuityModelingPrediction: FOR BEST-FITTING COLOR-OPPONENT MODEL
% acuityModelingPredictionLum: FOR BEST-FITTING LUMINANCE MODEL
fileStr = 'acuityModelingPrediction';

figure;
set(gcf,'Position',[176 273 1309 669]);
for j = 1:length(subjNumAll) % FOR EACH SUBJECT
    % GRAB SUBJECT NUMBER
    subjNum = subjNumAll(j);
    
    % LOAD PRE-GENERATED DATA
    load(fullfile(dataPath,'data','acuityModeling',[fileStr 'S' num2str(subjNum) '.mat']));

    scaleFac = 0.816; % ACCOUNTS FOR OPTOTUNE TO ACTUAL DEFOCUS CONVERSION
    
    % 'DEPTH-OF-FOCUS' PARAMETER (SET 0 NOW)
    shiftVals = 0;
    
    % SHIFTING PREDICTED D-PRIME CURVE AND REGRESSING AGAINST ACTUAL
    % D-PRIME VALUES TO GET THE BEST FIT
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
    dprimeScaleBestFit = dprimeScaleTmp(indMinShift);
    
    subplot(2,4,j);
    hold on;
    plot(shiftValBestFit+defocusForStim+modelPrediction875nmPurpleAt2pt5-2.5,normcdf(dprimeMetric.*dprimeScaleBestFit/2),'-','Color',[0.56 0 1],'LineWidth',1);
    errorbar(unqFocDst.*scaleFac,normcdf(dprime/2),(normcdf(dprime/2)-normcdf(dprimeCI(1,:)/2)),(normcdf(dprimeCI(2,:)/2)-normcdf(dprime/2)),'o','Color',[0.56 0 1],'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',8);
    xlabel('Relative Distance (D)');
    ylabel('Proportion Correct');
    set(gca,'FontSize',12);
    set(gca,'Box','on');
    set(gca,'XTick',[-1 -0.5 0 0.5 1]);
    axis square;
    text(-1,0.3,['S' num2str(j)],'FontSize',18);
    ylim([0.2 1.05]);
    xlim([-1.21 1.21]);
    
    % DETERMINE PREDICTED PEAK PERFORMANCE AFTER FITTING DEPTH-OF-FOCUS
    % PARAMETER
    stimDistanceSmp = 1.2:0.01:3.8;
    dprimeMetricSmooth = interp1(defocusForStim+modelPrediction875nmPurpleAt2pt5+shiftVals(indMinShift),dprimeMetric,stimDistanceSmp,'spline');
    [~,indPeak] = max(dprimeMetricSmooth);
    peakLocModelPrediction = stimDistanceSmp(indPeak);
    peakLocModelPredictionAll(j) = peakLocModelPrediction;
end
