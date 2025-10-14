%% GENERATE MODEL PREDICTIONS TOGETHER WITH ACTUAL DATA

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
% dataPath = 'G:\Shared drives\CIVO_BVAMS\';
subjNum = [1 3 5 10 16 17 18 20];
wvMeanAll = [];
wvPredAll = [];
aicAll = [];
dfPredPurpleAll = [];
wLMminAll = [];
wLpropMinAll = [];

for i = 1:length(subjNum)
    [aic, pFit, wvMean, wvPred, dfPredPurple, wLMmin, wLpropMin] = ARCtestWvInFocusMeanZspatFilterLMSeffectPlotWave(subjNum(i),'LM',dataPath);
    wvMeanAll(:,:,i) = wvMean;
    wvPredAll(:,:,i) = wvPred;
    aicAll(i) = aic;
    dfPredPurpleAll(i) = dfPredPurple;
    wLMminAll(i) = wLMmin;
    wLpropMinAll(i) = wLpropMin;
end 

%% MAKE FIGURE 4B IN THE PAPER

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];

load([foldername 'wvMeanAndPredLM.mat']);

symbDist = 'sod';
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
                         1.00 1.00 1.00];

defocusDiffDue2LCA = [];
predDiffDue2LCA = [];
subjNum = [1 3 5 10 16 17 18 20];

figure;
set(gcf,'Position',[180 156 959 574]);
subplot(1,2,1);
hold on;
for i = 1:size(wvMeanAll,2)
    wvMeanTmp = squeeze(wvMeanAll(:,i,:));
    wvPredTmp = squeeze(wvPredAll(:,i,:));
    plot(mean(wvPredTmp(1:5,:),2),'k-','LineWidth',1);
    for j = 1:5
        errorbar(j,mean(wvMeanTmp(j,:),2),1.96.*std(wvMeanTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
        for k = 1:size(wvMeanTmp,2)
            defocusDiffDue2LCA(i,j,k) = humanWaveDefocusARC(wvMeanTmp(1,k),wvMeanTmp(j,k),subjNum(k));
            predDiffDue2LCA(i,j,k) = humanWaveDefocusARC(wvPredTmp(1,k),wvPredTmp(j,k),subjNum(k));
        end
    end
end
axis square;
set(gca,'FontSize',15);
set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
set(gca,'Box','on');
xlim([0.5 5.5]);
ylim([450 680]);
xlabel('Red-blue ratio');
ylabel('Wavelength in focus (nm)');
subplot(1,2,2);
hold on;
for i = 1:size(wvMeanAll,2)
    wvMeanTmp = squeeze(wvMeanAll(:,i,:));
    wvPredTmp = squeeze(wvPredAll(:,i,:));
    plot(mean(wvPredTmp(6:10,:),2),'k-','LineWidth',1);
    for j = 6:10
        errorbar(j-5,mean(wvMeanTmp(j,:),2),1.96.*std(wvMeanTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
    end
end
axis square;
xlim([0.5 5.5]);
ylim([450 680]);
set(gca,'Box','on');
set(gca,'FontSize',15);
set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
xlabel('Red-blue ratio');

defocusDiffDue2LCArb = mean(squeeze(defocusDiffDue2LCA(:,5,:)));
predDiffDue2LCArb = mean(squeeze(predDiffDue2LCA(:,5,:)));

defocusDiffDue2LCArbSE = 1.96.*std(defocusDiffDue2LCArb)./sqrt(8);
predDiffDue2LCArbSE = 1.96.*std(predDiffDue2LCArb)./sqrt(8);
%%

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/wvMeanAndPredLM.mat');

symbDist = 'sod';
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
                         1.00 1.00 1.00];

figure;
set(gcf,'Position',[381 451 900 444])
subplot(1,2,1);
hold on;
for i = 1:size(wvMeanAll,2)
    wvMeanTmp = squeeze(wvMeanAll(:,i,:));
    wvPredTmp = squeeze(wvPredAll(:,i,:));
    plot(mean(wvMeanTmp(1:5,:),2),'k-','LineWidth',1);
    for j = 1:5
        errorbar(j,mean(wvMeanTmp(j,:),2),std(wvMeanTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
    end
end
% axis square;
set(gca,'FontSize',15);
set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
xlim([0.5 5.5]);
ylim([450 670]);
axis square;
box on;
xlabel('Red-blue ratio');
ylabel('Wavelength in focus (nm)');
subplot(1,2,2);
hold on;
for i = 1:size(wvMeanAll,2)
    wvMeanTmp = squeeze(wvMeanAll(:,i,:));
    wvPredTmp = squeeze(wvPredAll(:,i,:));
    plot(mean(wvMeanTmp(6:10,:),2),'k-','LineWidth',1);
    for j = 6:10
        errorbar(j-5,mean(wvMeanTmp(j,:),2),std(wvMeanTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
    end
end
% axis square;
xlim([0.5 5.5]);
ylim([450 670]);
axis square;
box on;
set(gca,'FontSize',15);
set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
xlabel('Red-blue ratio');
