%% MAKE FIGURE 5D

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

% plotSymb = '>^<vsdho'; % SYMBOLS FOR PLOTTING INDIVIDUALS
plotSymb = 'oooooooo'; % GIVING EACH INDIVIDUAL SAME SYMBOL

% acuityModelingPredictionLM: FOR BEST-FITTING LUMINANCE MODEL
[RMSELM, corrPredLM, pLM, peakLocModelPredictionAllLM, peakLocActualAllLM] = ARCnlz_Fig5Dhelper(dataPath,'acuityModelingPredictionLM');
% acuityModelingPrediction: FOR BEST-FITTING COLOR-OPPONENT MODEL
[RMSE, corrPred, p, peakLocModelPredictionAll, peakLocActualAll] = ARCnlz_Fig5Dhelper(dataPath,'acuityModelingPrediction');

% EXPLICITLY REPRESENT ERRORS OF EACH MODEL
errorAll = abs(peakLocModelPredictionAll-peakLocActualAll);
errorAllLM = abs(peakLocModelPredictionAllLM-peakLocActualAllLM);
% T-TEST ON ERRORS
[h,pError,~,stats] = ttest(errorAll-errorAllLM);

figure;
hold on;
for i = 1:length(peakLocModelPredictionAllLM)
   plot(peakLocModelPredictionAllLM(i),peakLocActualAllLM(i),['k' plotSymb(i)],'MarkerFaceColor',[1 1 1],'MarkerSize',15);
end
for i = 1:length(peakLocModelPredictionAll)
   plot(peakLocModelPredictionAll(i),peakLocActualAll(i),['k' plotSymb(i)],'MarkerFaceColor',[0.56 0 1],'MarkerSize',15);
end
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
legend({['\rho_c = ' num2str(corrPred,2)] ['\rho_l = ' num2str(corrPredLM,2)]},'Location','SouthEast');

% DISPLAY STATS

display('--------------------------------');
display(['Chromatic model correlation =  ' num2str(corrPred) ', p = ' num2str(p,3) newline]);
display(['Luminance model correlation = ' num2str(corrPredLM) ', p = ' num2str(pLM,3) newline]);
display(['Chromatic model RMSE =  ' num2str(RMSE,3) newline]);
display(['Chromatic model RMSE =  ' num2str(RMSELM,3) newline]);
display(['t(' num2str(stats.df) ') = ' num2str(stats.tstat,3) ', p = ' num2str(pError) newline]);