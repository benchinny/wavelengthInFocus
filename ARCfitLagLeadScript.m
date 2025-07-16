%%

subjNumAll = [1 3 5 10 16 17 18 20];
% subjNumAll = [10];
bSpatFilter = false;

SvaluesAll = [0 0 -1];
loadStr = {'0' '0' '-10'};    

if bSpatFilter
    coneWeightsFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsErrorSpatFilter/';
    figTag = 'filtered';
else
    coneWeightsFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsError/';
    figTag = 'unfiltered';
end

rmsAllSubj = [];
rmsBalanceSubj = [];
rmsMeanAllSubj = [];
rmsMeanBalanceSubj = [];
pFitMeanBalanceAll = [];
pFitMeanAllSubj = [];

for i = 1:length(subjNumAll)
    subjNum = subjNumAll(i);
    if subjNum==10
        RMSEall = zeros([2 2 2]);
        blockNums = 3:8;
        ind2examine = 1:3;     
        coordinates2examine = [1 1; 1 2; 1 1];
    elseif subjNum==3
        RMSEall = zeros([2 2 2]);
        blockNums = 12:17;
        ind2examine = 1:3;    
        coordinates2examine = [1 1; 1 2; 1 1];
    elseif subjNum==1
        RMSEall = zeros([2 2 2]);
        ind2examine = 1:3;
        blockNums = 11:16;
        coordinates2examine = [1 1; 1 2; 1 1];
    elseif subjNum==5
        RMSEall = zeros([2 2 2]);
        ind2examine = 1:3;
        blockNums = 3:8;
        coordinates2examine = [1 1; 1 2; 1 1];
    elseif subjNum==9
        RMSEall = zeros([2 2 2]);
        ind2examine = 1:3;
        blockNums = 2:7;
        coordinates2examine = [1 1; 1 2; 1 1];
    elseif subjNum==16
        RMSEall = zeros([2 2 2]);
        ind2examine = 1:3;
        blockNums = 2:7;
        coordinates2examine = [1 1; 1 2; 1 1];
    elseif subjNum==17
        RMSEall = zeros([2 2 2]);
        ind2examine = 1:3;
        blockNums = 2:7;
        coordinates2examine = [1 1; 1 2; 1 1];
    elseif subjNum==18
        RMSEall = zeros([2 2 2]);
        ind2examine = 1:3;
        blockNums = 2:7;
        coordinates2examine = [1 1; 1 2; 1 1];
    elseif subjNum==20    
        RMSEall = zeros([2 2 2]);
        ind2examine = 1:3;
        blockNums = 2:7;
        coordinates2examine = [1 1; 1 2; 1 1];
    end
    rgbAll = [];
    optDistStim = [];
    
    for j = 1:length(blockNums)
        AFCp = ARCloadFileBVAMS(subjNum+10,blockNums(j));
        rgbAll = [rgbAll; AFCp.rgb100];
        optDistStim = [optDistStim; AFCp.meanv00./1.2255];
    end
    
    rgbLumNorm = [];
    rgbLumNorm(:,1) = (rgbAll(:,1).^2.5)./0.2442;
    rgbLumNorm(:,2) = (rgbAll(:,2).^2.7)./0.1037;
    rgbLumNorm(:,3) = (rgbAll(:,3).^2.3)./1;
    rgbLumNorm(rgbLumNorm>1) = 1;

    for j = ind2examine
        load([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResults' num2str(loadStr{j}) '.mat']);
        RMSEall(:,:,j) = RMSE(:,:,end);
    end

    for j = ind2examine
        RMSEtmp = squeeze(RMSEall(:,:,j)); 
        load([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResults' num2str(loadStr{j}) '.mat']);
        predD = squeeze(defocus875predAll(:,coordinates2examine(j,1),coordinates2examine(j,2),end))';
        predD = predD(:);
        actD = squeeze(defocus875all(:,coordinates2examine(j,1),coordinates2examine(j,2),end))';
        actD = actD(:); 
        [pFit,rms] = ARCfitLagLead(predD,actD,optDistStim);
        pFitAll(:,j) = pFit';
        predDall(:,j) = predD;
        actDall(:,j) = actD; 
        rmsAll(:,j) = rms;
        rmsAllSubj(i,j) = rms;
    end
    
    [pFitBalance,rmsBalance] = ARCfitLagLead(optDistStim,actD,optDistStim);
    rmsBalanceSubj(i,:) = rmsBalance;

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
    
    pFitMeanAll = [];
    rmsMeanAll = [];
    defocusAt875meanPredAll = [];
    defocusAt875meanAll = [];
    defocusAt875meanPredBalanceAll = [];

    figure;
    set(gcf,'Position',[95 155 1089 751]);
    for j = 1:length(optDistToCheckAll)
        subplot(2,2,j);
        hold on;
        optDistToCheck = optDistToCheckAll(j);
        for k = 1:size(conditionsOrderedNorm,1)
            ind = abs(rgbLumNorm(:,1)-conditionsOrderedNorm(k,1))<0.01 & ...
                  abs(rgbLumNorm(:,2)-conditionsOrderedNorm(k,2))<0.01 & ...
                  abs(rgbLumNorm(:,3)-conditionsOrderedNorm(k,3))<0.01 & ...
                  abs(optDistStim-optDistToCheck)<0.01 & ...
                  indGoodOpt;
            if k<size(conditionsOrderedNorm,1)
                plot(k.*ones([sum(ind) 1]),actD(ind),'o','Color',conditionsOrderedNorm(k,:),'MarkerFaceColor',conditionsOrderedNorm(k,:));
            else
                plot(k.*ones([sum(ind) 1]),actD(ind),'o','Color','k','MarkerFaceColor','k');
            end
            defocusAt875mean(k) = mean(actD(ind));
            defocusAt875meanPred(k,1) = mean(predDall(ind,1));
            defocusAt875meanPred(k,2) = mean(predDall(ind,2));
            defocusAt875meanPred(k,3) = mean(predDall(ind,3));
        end

        for k = 1:size(predDall,2)
            [pFitMean,rmsMean] = ARCfitLagLead(defocusAt875meanPred(:,k),defocusAt875mean',optDistToCheck.*ones(size(defocusAt875mean')));
            pFitMeanAll(j,k) = pFitMean;
        end

        [pFitMeanBalance(j,:),rmsMeanBalance(j,:)] = ARCfitLagLead(optDistToCheck.*ones(size(defocusAt875mean')),defocusAt875mean',optDistToCheck.*ones(size(defocusAt875mean')));

        defocusAt875meanAll = [defocusAt875meanAll; defocusAt875mean'];
        defocusAt875meanPredAll = [defocusAt875meanPredAll; ...
                                   defocusAt875meanPred(:,1)-pFitMeanAll(j,1) ...
                                   defocusAt875meanPred(:,2)-pFitMeanAll(j,2) ...
                                   defocusAt875meanPred(:,3)-pFitMeanAll(j,3)];
        defocusAt875meanPredBalanceAll = [defocusAt875meanPredBalanceAll; optDistToCheck.*ones(size(defocusAt875mean'))-pFitMeanBalance(j)];

        plot(defocusAt875mean(1:5),'-','LineWidth',1,'Color',0.7.*[1 1 1]);
        plot(6:11,defocusAt875mean(6:11),'-','LineWidth',1,'Color',0.7*[1 1 1]);

        plot(defocusAt875meanPred(1:5,1)-pFitMeanAll(j,1),'-','Color',0.*[1 1 1],'LineWidth',1.5);
        plot(6:11,defocusAt875meanPred(6:11,1)-pFitMeanAll(j,1),'-','Color',0.*[1 1 1],'LineWidth',1.5);
        plot(defocusAt875meanPred(1:5,2)-pFitMeanAll(j,2),'-','Color',1.*[0 1 0],'LineWidth',1.5);
        plot(6:11,defocusAt875meanPred(6:11,2)-pFitMeanAll(j,2),'-','Color',1.*[0 1 0],'LineWidth',1.5);
        plot(defocusAt875meanPred(1:5,3)-pFitMeanAll(j,3),'-','Color',1.*[0 0 1],'LineWidth',1.5);
        plot(6:11,defocusAt875meanPred(6:11,3)-pFitMeanAll(j,3),'-','Color',1.*[0 0 1],'LineWidth',1.5); 

        % [defocusAt875meanPred(:,1)-pFitMeanAll(j,1) defocusAt875mean']

        plot([0 11],defocusAt875mean(11).*[1 1],'k--','LineWidth',1);
        xlim([0 11]);
        % ylim(mean(actD(indDist))+[-0.6 0.6]);
        title(['Subject ' num2str(subjNum) ', Optical Distances = ' num2str(optDistToCheck)]);
        plot(5.5.*[1 1],ylim,'k-');
        set(gca,'FontSize',15);
        xlabel('Condition');
        ylabel('Defocus at 875nm');
    end
    rmsMeanAllSubj(i,1) = sqrt(mean((defocusAt875meanAll-defocusAt875meanPredAll(:,1)).^2));
    rmsMeanAllSubj(i,2) = sqrt(mean((defocusAt875meanAll-defocusAt875meanPredAll(:,2)).^2));
    rmsMeanAllSubj(i,3) = sqrt(mean((defocusAt875meanAll-defocusAt875meanPredAll(:,3)).^2));
    rmsMeanBalanceSubj(i,:) = sqrt(mean((defocusAt875meanAll-defocusAt875meanPredBalanceAll).^2));
    pFitMeanBalanceAll(:,i) = pFitMeanBalance;
    pFitMeanAllSubj(:,:,i) = pFitMeanAll;
    % saveas(gcf,[coneWeightsFolder 'colorMechPredictionsS' num2str(subjNum) figTag],'png');
end
