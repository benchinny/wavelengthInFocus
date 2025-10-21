%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
subjNumAll = [11 13 15 20 26 27 28 30]-10;
aicAll = [];

for i = 1:length(subjNumAll)
    [aic, pFit, wvMean, wvPred] = ARCtestWvInFocusMeanZdeltaPassPlotStack(subjNumAll(i),'LM',4,dataPath);
    aicAll(i) = aic;
    wvMeanAll(:,:,i) = wvMean;
    wvPredAll(:,:,i) = wvPred;
    pFitAll(:,i) = pFit;    
end

%%

if ispc
slash = '\';
else
slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];
save([foldername 'wvMeanAndPredDeltaPass2.mat'],'wvMeanAll','wvPredAll','aicAll','pFitAll');

%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];

load([foldername 'wvMeanAndPredDeltaPass2.mat']);

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

aicAll2 = [-24.2269 -33.3423 -4.7110 -14.0518 -42.0595 -44.0608 -15.9615 1.4181];

