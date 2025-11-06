%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];
plotSymAll = 'xosd>^<v';

load([foldername 'wvMeanAndPredDonutx2.mat']);
dfPredPurpleAllBY = dfPredPurpleAll;
aicBY = aicAll;
load([foldername 'wvMeanAndPredLminusM.mat']);
dfPredPurpleAllRG = dfPredPurpleAll;
aicRG = aicAll;
load([foldername 'wvMeanAndPredLM.mat']);
dfPredPurpleAllLum = dfPredPurpleAll;
aicLum = aicAll;

figure; 
set(gcf,'Position',[1 1 1132 539]);
subplot(1,3,1);
hold on;
for i = 1:length(dfPredPurpleAllBY)
   plot(-dfPredPurpleAllBY(i),-dfPredPurpleAllRG(i),['k' plotSymAll(i)],'MarkerSize',10, ...
       'MarkerFaceColor','w'); 
end
axis square; 
formatFigure('Blue-Yellow','Red-Green',['\rho = ' num2str(corr(dfPredPurpleAllBY',dfPredPurpleAllRG'),2)]);
xlim([1 2]);
ylim([1 2]);
subplot(1,3,2);
hold on;
for i = 1:length(dfPredPurpleAllLum)
   plot(-dfPredPurpleAllLum(i),-dfPredPurpleAllBY(i),['k' plotSymAll(i)],'MarkerSize',10, ...
       'MarkerFaceColor','w'); 
end
axis square; 
formatFigure('Luminance','Blue-yellow',['\rho = ' num2str(corr(dfPredPurpleAllLum',dfPredPurpleAllBY'),2)]);
xlim([1 2]);
ylim([1 2]);
subplot(1,3,3);
hold on;
for i = 1:length(dfPredPurpleAllLum)
   plot(-dfPredPurpleAllLum(i),-dfPredPurpleAllRG(i),['k' plotSymAll(i)],'MarkerSize',10, ...
       'MarkerFaceColor','w'); 
end
axis square; 
formatFigure('Luminance','Red-green',['\rho = ' num2str(corr(dfPredPurpleAllLum',dfPredPurpleAllRG'),2)]);
xlim([1 2]);
ylim([1 2]);

%% COMPARISON AFTER DEPTH-OF-FOCUS PARAMETER IS APPLIED

plotSymAll = 'xosd>^<v';
dfFitPurpleAllBY = [1.81 2.06 2.34 2.20 2.18 1.96 1.81 2.09];
dfFitPurpleAllRG = [1.81 2.06 2.35 2.17 2.22 1.97 1.83 2.10];
dfFitPurpleAllLum = [1.89 2.06 2.40 2.19 1.52 1.97 1.81 2.06];

figure; 
set(gcf,'Position',[1 1 1132 539]);
subplot(1,3,1);
hold on;
for i = 1:length(dfFitPurpleAllBY)
    plot(dfFitPurpleAllRG(i),dfFitPurpleAllBY(i),['k' plotSymAll(i)],'MarkerSize',10,'MarkerFaceColor','w');
end
axis square;
formatFigure('Red-green','Blue-yellow');
xlim([1.25 2.75]);
ylim([1.25 2.75]);
subplot(1,3,2);
hold on;
for i = 1:length(dfFitPurpleAllLum)
    plot(dfFitPurpleAllLum(i),dfFitPurpleAllBY(i),['k' plotSymAll(i)],'MarkerSize',10,'MarkerFaceColor','w');
end
axis square;
formatFigure('Luminance','Blue-yellow','After depth-of-focus parameter');
xlim([1.25 2.75]);
ylim([1.25 2.75]);
subplot(1,3,3);
hold on;
for i = 1:length(dfFitPurpleAllLum)
    plot(dfFitPurpleAllLum(i),dfFitPurpleAllRG(i),['k' plotSymAll(i)],'MarkerSize',10,'MarkerFaceColor','w');
end
axis square;
formatFigure('Luminance','Red-green');
xlim([1.25 2.75]);
ylim([1.25 2.75]);

%% FIGURE 5C WITH SYMBOLS

plotSymAll = 'xosd>^<v';
dprimeRatioAll = [2.41 1.000 1.000 1.41 1.000 1.98 1.76 1.023];

figure;
% boxplot(dprimeRatioAll,'LineWidth',1,'FaceColor',[0.56 0 1]);
boxplot(dprimeRatioAll);
set(gca,'FontSize',15);
hold on;
for i = 1:length(dprimeRatioAll)
    if abs(dprimeRatioAll(i)-1)<0.01
        plot(1+i*0.1,dprimeRatioAll(i),plotSymAll(i),'MarkerSize',10,'MarkerFaceColor',[0.56 0 1],'Color',[0.56 0 1]);
    else
        plot(1,dprimeRatioAll(i),plotSymAll(i),'MarkerSize',10,'MarkerFaceColor',[0.56 0 1],'Color',[0.56 0 1]);
    end
end
% set(gca,'XTick',1:8);
% set(gca,'XTickLabel',{'S1' 'S3' 'S5' 'S10' 'S16' 'S17' 'S18' 'S20'});
% xlabel('Subject');
ylabel('d''_{max}/d''_{0}');
% ylim([0 1]);