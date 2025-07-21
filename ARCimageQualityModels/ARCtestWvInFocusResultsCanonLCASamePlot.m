%%

clear;

%%

subjNum = 30;

% 1: L+M
% 2: L-M
% 3: L+M-S
% 4: All
% 5: Best fit

SvaluesAll = [0 0 -1];
loadStr = {'0' '0' '-10'};    
optDistNames = {'1pt5' '2pt5' '3pt5'};

if subjNum==20
    % RMSEall = zeros([11 11 11]);
    % SvaluesAll = [-1 -0.8 -0.6 -0.4 -0.2 0.0 0.2 0.4 0.6 0.8 1];
    % loadStr = {'-10' '-8' '-6' '-4' '-2' '0' '2' '4' '6' '8' '10'};
    % blockNums = 3:8;
    % ind2examine = 1:11;

    % RMSEall = zeros([11 11 1]);
    % SvaluesAll = [-0.8];
    % loadStr = {'-8'};
    % blockNums = 3:8;
    % ind2examine = 1;    
    RMSEall = zeros([2 2 2]);
    blockNums = 3:8;
    ind2examine = 1:3;     
    coordinates2examine = [1 1; 1 2; 1 1];
elseif subjNum==13
    RMSEall = zeros([2 2 2]);
    blockNums = 12:17;
    ind2examine = 1:3;    
    coordinates2examine = [1 1; 1 2; 1 1];
elseif subjNum==11
    RMSEall = zeros([2 2 2]);
    ind2examine = 1:3;
    blockNums = 11:16;
    coordinates2examine = [1 1; 1 2; 1 1];
elseif subjNum==15
    RMSEall = zeros([2 2 2]);
    ind2examine = 1:3;
    blockNums = 3:8;
    coordinates2examine = [1 1; 1 2; 1 1];
elseif subjNum==19
    RMSEall = zeros([2 2 2]);
    ind2examine = 1:3;
    blockNums = 2:7;
    coordinates2examine = [1 1; 1 2; 1 1];
elseif subjNum==26
    RMSEall = zeros([2 2 2]);
    ind2examine = 1:3;
    blockNums = 2:7;
    coordinates2examine = [1 1; 1 2; 1 1];
elseif subjNum==27
    RMSEall = zeros([2 2 2]);
    ind2examine = 1:3;
    blockNums = 2:7;
    coordinates2examine = [1 1; 1 2; 1 1];
elseif subjNum==28
    RMSEall = zeros([2 2 2]);
    ind2examine = 1:3;
    blockNums = 2:7;
    coordinates2examine = [1 1; 1 2; 1 1];
elseif subjNum==30    
    RMSEall = zeros([2 2 2]);
    ind2examine = 1:3;
    blockNums = 2:7;
    coordinates2examine = [1 1; 1 2; 1 1];
end

folderPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsErrorOrigLCA/';
rgbAll = [];
optDistStim = [];

for i = 1:length(blockNums)
    AFCp = ARCloadFileBVAMS(subjNum,blockNums(i));
    rgbAll = [rgbAll; AFCp.rgb100];
    optDistStim = [optDistStim; AFCp.meanv00./1.2255];
end

rgbLumNorm = [];
rgbLumNorm(:,1) = (rgbAll(:,1).^2.5)./0.2442;
rgbLumNorm(:,2) = (rgbAll(:,2).^2.7)./0.1037;
rgbLumNorm(:,3) = (rgbAll(:,3).^2.3)./1;
rgbLumNorm(rgbLumNorm>1) = 1;

for i = ind2examine
    load([folderPath 'S' num2str(subjNum-10) 'wvInFocusModelResults' num2str(loadStr{i}) '.mat']);
    RMSEall(:,:,i) = RMSE(:,:,end);
end

globalMinRMSE = min(RMSEall(:));

predDall = [];
actDall = [];

for i = ind2examine
    RMSEtmp = squeeze(RMSEall(:,:,i)); 

    figure; 
    imagesc(RMSEtmp);
    clim([min(RMSEall(:)) max(RMSEall(:))]);
    xlabel('M weight');
    ylabel('L weight');
    axis square;
    set(gca,'FontSize',12);
    set(gca,'XTick',1:11);
    set(gca,'XTickLabel',{'-1' '-0.8' '-0.6' '-0.4' '-0.2' '0' '0.2' '0.4' '0.6' '0.8' '1'});
    set(gca,'YTickLabel',{'-1' '-0.8' '-0.6' '-0.4' '-0.2' '0' '0.2' '0.4' '0.6' '0.8' '1'});
    title(['S weight = ' num2str(SvaluesAll(i))]);   
    minRMSE = min(RMSEtmp(:));
    if minRMSE==globalMinRMSE
        minCoordinates = RMSEtmp == minRMSE;
        SvalueMin = SvaluesAll(i);
    end
    
    load([folderPath 'S' num2str(subjNum-10) 'wvInFocusModelResults' num2str(loadStr{i}) '.mat']);
    predD = squeeze(defocus875predAll(:,coordinates2examine(i,1),coordinates2examine(i,2),end))';
    predD = predD(:);
    actD = squeeze(defocus875all(:,coordinates2examine(i,1),coordinates2examine(i,2),end))';
    actD = actD(:);
    figure; 
    hold on;
    for j = 1:length(predD)
        plot(predD(j),actD(j),'.','MarkerSize',12,'Color',rgbLumNorm(j,:)); 
    end
    plot([0 3.5],[0 3.5],'k--','LineWidth',1);
    xlim([0 3.5]);
    ylim([0 3.5]);
    axis square;
    formatFigure('Predicted Defocus (D)','Actual Defocus (D)');    

    predDall(:,i) = predD;
    actDall(:,i) = actD;
end

%% PLOTTING ALL TRIAL MEANS PER CONDITION AND DISTANCE

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

figPositions = [14 493 560 420; ...
                544 496 560 420; ...
                1079 498 560 420; ...
                ];
optDistToCheckAll = [1.5 2.5 3.5];
indGoodOpt = abs(actD+1-optDistStim)<1;

figure;
set(gcf,'Position',[95 155 1089 751]);
for j = 1:length(optDistToCheckAll)
    subplot(2,2,j);
    hold on;
    optDistToCheck = optDistToCheckAll(j);
    indDist = abs(optDistStim-optDistToCheck)<0.01;
    for i = 1:size(conditionsOrderedNorm,1)
        ind = abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
              abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
              abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01 & ...
              abs(optDistStim-optDistToCheck)<0.01 & ...
              indGoodOpt;
        if i<size(conditionsOrderedNorm,1)
            plot(i.*ones([sum(ind) 1]),actD(ind),'o','Color',conditionsOrderedNorm(i,:),'MarkerFaceColor',conditionsOrderedNorm(i,:));
        else
            plot(i.*ones([sum(ind) 1]),actD(ind),'o','Color','k','MarkerFaceColor','k');
        end
        defocusAt875mean(i) = mean(actD(ind));
        defocusAt875meanPred(i,1) = mean(predDall(ind,1));
        defocusAt875meanPred(i,2) = mean(predDall(ind,2));
        defocusAt875meanPred(i,3) = mean(predDall(ind,3));
    end
    plot(defocusAt875mean(1:5),'k-');
    plot(6:10,defocusAt875mean(6:10),'k-');
    plot(defocusAt875meanPred(1:5,1),'-','Color',0.*[1 1 1],'LineWidth',1.5);
    plot(6:10,defocusAt875meanPred(6:10,1),'-','Color',0.*[1 1 1],'LineWidth',1.5);
    plot(defocusAt875meanPred(1:5,2),'-','Color',1.*[0 1 0],'LineWidth',1.5);
    plot(6:10,defocusAt875meanPred(6:10,2),'-','Color',1.*[0 1 0],'LineWidth',1.5);
    plot(defocusAt875meanPred(1:5,3),'-','Color',1.*[0 0 1],'LineWidth',1.5);
    plot(6:10,defocusAt875meanPred(6:10,3),'-','Color',1.*[0 0 1],'LineWidth',1.5);    
    plot([0 11],defocusAt875mean(11).*[1 1],'k--','LineWidth',1);
    xlim([0 11]);
    % ylim(mean(actD(indDist))+[-0.6 0.6]);
    title(['Subject ' num2str(subjNum-10) ', Optical Distances = ' num2str(optDistToCheck)]);
    plot(5.5.*[1 1],ylim,'k-');
    set(gca,'FontSize',15);
    xlabel('Condition');
    ylabel('Defocus at 875nm');
end
saveas(gcf,['/Users/benjaminchin/Documents/ARchromaScraps/colorMechPredictionsS' num2str(subjNum-10) 'mechAllOrigLCA'],'png');

% %% PLOTTING ALL TRIAL MEANS PER CONDITION AND DISTANCE
% 
% conditionsOrderedNorm = [0.25 0.00 1.00; ...
%                          0.50 0.00 1.00; ...
%                          1.00 0.00 1.00; ...
%                          1.00 0.00 0.50; ...
%                          1.00 0.00 0.25; ...
%                          0.25 0.50 1.00; ...
%                          0.50 0.50 1.00; ...
%                          1.00 0.50 1.00; ...
%                          1.00 0.50 0.50; ...
%                          1.00 0.50 0.25; ...
%                          1.00 1.00 1.00];
% 
% figPositions = [15 78 560 420; ...
%                 553 81 560 420; ...
%                 1079 85 560 420; ...
%                 ];
% optDistToCheckAll = [1.5 2.5 3.5];
% indGoodOpt = abs(actD+1-optDistStim)<1;
% 
% for j = 1:length(optDistToCheckAll)
%     figure;
%     set(gcf,'Position',figPositions(j,:));
%     hold on;
%     optDistToCheck = optDistToCheckAll(j);
%     indDist = abs(optDistStim-optDistToCheck)<0.01;
%     for i = 1:size(conditionsOrderedNorm,1)
%         ind = abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
%               abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
%               abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01 & ...
%               abs(optDistStim-optDistToCheck)<0.01 & ...
%               indGoodOpt;
%         if i<size(conditionsOrderedNorm,1)
%             plot(i.*ones([sum(ind) 1]),predD(ind),'o','Color',conditionsOrderedNorm(i,:),'MarkerFaceColor',conditionsOrderedNorm(i,:));
%         else
%             plot(i.*ones([sum(ind) 1]),predD(ind),'o','Color','k','MarkerFaceColor','k');
%         end
%         defocusAt875mean(i) = mean(predD(ind));
%     end
%     plot(defocusAt875mean(1:5),'k-');
%     plot(6:10,defocusAt875mean(6:10),'k-');
%     plot([0 11],defocusAt875mean(11).*[1 1],'k--','LineWidth',1);
%     xlim([0 11]);
%     % ylim(mean(actD(indDist))+[-0.6 0.6]);
%     title(['Subject ' num2str(subjNum-10) ', Optical Distances = ' num2str(optDistToCheck)]);
%     plot(5.5.*[1 1],ylim,'k-');
%     set(gca,'FontSize',15);
%     xlabel('Condition');
%     ylabel('Defocus at 875nm');
% end
