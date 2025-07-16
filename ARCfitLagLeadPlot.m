%%

clear all;

%%

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsError/rmsMeanAllSubjComb.mat');

%%

subjNumAll = [1 3 5 10 16 17 18 20];

figure;
set(gcf,'Position',[123 316 1476 621]);
for i = 1:size(rmsMeanAllSubjComb,1)
    subplot(2,4,i);
    hold on;
    bar(1,rmsMeanAllSubjComb(i,1),'FaceColor','w','LineWidth',1);
    bar(2,rmsMeanAllSubjComb(i,2),'FaceColor',[0 1 0],'LineWidth',1);
    bar(3,rmsMeanAllSubjComb(i,3),'FaceColor',[0 0 1],'LineWidth',1);
    bar(4,rmsMeanAllSubjComb(i,4),'FaceColor',0.5.*[1 1 1],'LineWidth',1);
    bar(5,rmsMeanAllSubjComb(i,5),'FaceColor',[0 0.5 0],'LineWidth',1);
    bar(6,rmsMeanAllSubjComb(i,6),'FaceColor',[0 0 0.5],'LineWidth',1);   
    bar(7,rmsMeanAllSubjComb(i,7),'FaceColor',[0 0 0],'LineWidth',1);   
    set(gca,'FontSize',10);
    set(gca,'XTick',[1 2 3 7]);
    set(gca,'XTickLabel',{'L+M' 'L-M' 'L+M-S' 'Bal'});
    xlabel('Model');
    ylabel('RMS');
    title(['Subject ' num2str(subjNumAll(i))]);
end

%% PLOT LAGS AND LEADS

subjNumAll = [1 3 5 10 16 17 18 20];

figure;
set(gcf,'Position',[141 322 1398 658]);
for i = 1:size(pFitMeanAllSubj,3)
    subplot(2,4,i);
    hold on;
    plot([1.5 2.5 3.5],pFitMeanAllSubj(:,1,i),'-ko','LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',10);
    plot([1.5 2.5 3.5],pFitMeanAllSubj(:,2,i),'-go','LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',10);
    plot([1.5 2.5 3.5],pFitMeanAllSubj(:,3,i),'-bo','LineWidth',1.5,'MarkerFaceColor','w','MarkerSize',10);
    axis square;
    plot([1 4],[0 0],'k--','LineWidth',1);
    xlim([1 4]);
    set(gca,'FontSize',15);
    if i==5
        xlabel('Stimulus distance (D)');
    end
    if i==1
        ylabel('Prediction correction (D)');
    end
end

