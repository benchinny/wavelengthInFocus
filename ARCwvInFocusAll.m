%%

subjNum = [1 3 5 10 16 17 18 20];
wvMeanAll = [];
wvPredAll = [];

for i = 1:length(subjNum)
    [aic, pFit, wvMean, wvPred] = ARCtestWvInFocusMeanZspatFilterLMSeffectPlotWave(subjNum(i),'LminusM');
    wvMeanAll(:,:,i) = wvMean;
    wvPredAll(:,:,i) = wvPred;
end 

%%

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/wvMeanAndPred.mat');

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
set(gcf,'Position',[737 404 593 333]);
subplot(1,2,1);
hold on;
for i = 1:size(wvMeanAll,2)
    wvMeanTmp = squeeze(wvMeanAll(:,i,:));
    wvPredTmp = squeeze(wvPredAll(:,i,:));
    plot(mean(wvPredTmp(1:5,:),2),'k-','LineWidth',1);
    for j = 1:5
        errorbar(j,mean(wvMeanTmp(j,:),2),std(wvPredTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
    end
end
% axis square;
set(gca,'FontSize',15);
set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
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
        errorbar(j-5,mean(wvMeanTmp(j,:),2),std(wvPredTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
    end
end
% axis square;
xlim([0.5 5.5]);
ylim([450 680]);
set(gca,'FontSize',15);
set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
xlabel('Red-blue ratio');

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
        errorbar(j,mean(wvMeanTmp(j,:),2),std(wvPredTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
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
        errorbar(j-5,mean(wvMeanTmp(j,:),2),std(wvPredTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
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
