%% COMPARE WEIGHTING MODEL AND SWITCHING MODEL

sn = [18 23 24 26 27 28 29 30 32 33];

rhoSwitchAll = [];
rhoFullAll = [];
rhoNoColorAll = [];
weightsRBSall = [];
rhoColorAll = [];
rhoColorSwitchAll = [];
dAll = [];
aicSwitchAll = [];
aicLinAll = [];
aicNoColorAll = [];
weightsRBSciAll = [];

% for i = 1:length(sn)
%     [d, wS, rbThreshold, rhoSwitch, rhoColorSwitch, aicSwitch] = ARCnlz_linearSwitching(sn(i),false);
%     rhoSwitchAll(i) = rhoSwitch;
%     rhoColorSwitchAll(i) = rhoColorSwitch;
%     dAll(i) = d;
%     aicSwitchAll(i) = aicSwitch;
%     display(['Done with subject ' num2str(i)]);
% end

for i = 1:length(sn)
    [weightsRBS1_x, weightsRBS1_y, rhoFull, rhoNoColor, rhoColor, aicLin, aicNoColor, weightsRBSci, trialMeans, meanChange, deltaR, deltaB, deltaS] = ARCnlz_linearModelnobias(sn(i),false,1000,0);
    rhoFullAll(i) = rhoFull;
    rhoNoColorAll(i) = rhoNoColor;
    weightsRBSall(i,:) = weightsRBS1_x;
    rhoColorAll(i) = rhoColor;
    aicLinAll(i) = aicLin;
    aicNoColorAll(i) = aicNoColor;
    weightsRBSciAll(:,:,i) = weightsRBSci;
    meanChangeAll(:,i) = meanChange;
    deltaRall(:,i) = deltaR;
    deltaBall(:,i) = deltaB;
    deltaSall(:,i) = deltaS;
    display(['Done with subject ' num2str(i)]);
end

%%

maxLumCdm2 = 0.40;
RBratioAll = [];
figure;
set(gcf,'Position',[119 203 1437 718]);
for i = 1:length(sn)
    subplot(2,5,i);
%    plot([deltaR deltaB deltaS]*weightsRBS1_x,meanChangeXvec,'ko','LineWidth',1);
    deltaR = deltaRall(:,i);
    deltaB = deltaBall(:,i);
    deltaS = deltaSall(:,i);
    meanChange = meanChangeAll(:,i);
    weightsRBS1_x = weightsRBSall(i,:);
    predictionTmp = [deltaR deltaB deltaS]*weightsRBS1_x';
    hold on;
    for j = 1:length(predictionTmp)
        deltaRtmp = (1/maxLumCdm2)*deltaR(j);
        deltaBtmp = (1/maxLumCdm2)*deltaB(j);
        RBratio = 0.3.*(deltaRtmp-deltaBtmp+2);
        if RBratio<0
            RBratio = 0;
        end
        if RBratio>1
            RBratio = 1;
        end
        RBratioAll(j)=RBratio;
        plot(predictionTmp(j),meanChange(j),'o','LineWidth',1,'Color',[RBratio 0 1-RBratio],'MarkerFaceColor',[RBratio 0 1-RBratio]);
        % plot(predictionTmp(j),meanChangeXvec(j),'o','LineWidth',1,'Color','k','MarkerFaceColor','w');
    end    
    plot([-3 3],[-3 3],'k--');
    xlim([-3 3]);
    ylim([-3 3]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaA');
    ylabel('Measured \DeltaA');
    title(['Correlation = ' num2str(rhoFullAll(i),3)]);
    axis square;
end

%%

maxLumCdm2 = 0.40;
RBratioAll = [];
figure;
set(gcf,'Position',[119 203 1437 718]);
for i = 1:length(sn)
    subplot(2,5,i);
%    plot([deltaR deltaB deltaS]*weightsRBS1_x,meanChangeXvec,'ko','LineWidth',1);
    deltaR = deltaRall(:,i);
    deltaB = deltaBall(:,i);
    deltaS = deltaSall(:,i);
    meanChange = meanChangeAll(:,i);
    weightsRBS1_x = weightsRBSall(i,:);
    predictionTmp = [deltaS]*weightsRBS1_x(3)';
    hold on;
    for j = 1:length(predictionTmp)
        deltaRtmp = (1/maxLumCdm2)*deltaR(j);
        deltaBtmp = (1/maxLumCdm2)*deltaB(j);
        RBratio = 0.3.*(deltaRtmp-deltaBtmp+2);
        if RBratio<0
            RBratio = 0;
        end
        if RBratio>1
            RBratio = 1;
        end
        RBratioAll(j)=RBratio;
        plot(predictionTmp(j),meanChange(j),'o','LineWidth',1,'Color',[RBratio 0 1-RBratio],'MarkerFaceColor',[RBratio 0 1-RBratio]);
        % plot(predictionTmp(j),meanChangeXvec(j),'o','LineWidth',1,'Color','k','MarkerFaceColor','w');
    end    
    plot([-3 3],[-3 3],'k--');
    xlim([-3 3]);
    ylim([-3 3]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaA');
    ylabel('Measured \DeltaA');
    title(['Correlation = ' num2str(rhoNoColorAll(i),3)]);
    axis square;
end

%% COMPARE WEIGHTING MODEL AND SWITCHING MODEL

sn = [18 23 24 26 27 28 29 30 32 33];
% meanDiffPureRedVsBlue = [0.9607 1.0299 0.9122 0.1594 0.4769 0.5201 0.2107 0.2619 0.8158 0.3803];
meanDiffPureRedVsBlue = 0.8.*ones([1 10]);

rhoSwitchAll = [];
dAll = [];
aicSwitchAll = [];
wSall = [];
rbThresholdAll = [];
meanChangeAll = [];
deltaRB1all = [];
deltaRB2all = [];
deltaSall = [];

for i = 1:length(sn)
    [d, wS, rbThreshold, rhoSwitch, ~, aicSwitch,trialMeans,meanChange, deltaRB1, deltaRB2, deltaS] = ARCnlz_linearSwitching(sn(i),false,meanDiffPureRedVsBlue(i));
    rhoSwitchAll(i) = rhoSwitch;
    dAll(i) = d;
    wSall(i) = wS;
    rbThresholdAll(i) = rbThreshold;
    aicSwitchAll(i) = aicSwitch;
    meanChangeAll(:,i) = meanChange;
    deltaRB1all(:,i) = deltaRB1;
    deltaRB2all(:,i) = deltaRB2;
    deltaSall(:,i) = deltaS;    
    display(['Done with subject ' num2str(i)]);
end

%%

h = figure;
set(gcf,'Position',[119 203 1437 718]);
for i = 1:length(sn)
    subplot(2,5,i);
%    plot([deltaR deltaB deltaS]*weightsRBS1_x,meanChangeXvec,'ko','LineWidth',1);
    deltaRB1 = deltaRB1all(:,i);
    deltaRB2 = deltaRB2all(:,i);
    wS = wSall(i);
    d = dAll(i);
    rbThreshold = rbThresholdAll(i);
    deltaS = deltaSall(:,i);
    meanChange = meanChangeAll(:,i);
    c = zeros(size(deltaS));
    wS = wS.*ones(size(deltaS));
    
    c(deltaRB1<rbThreshold & deltaRB2>rbThreshold)=d;
    c(deltaRB1>rbThreshold & deltaRB2<rbThreshold)=-d;
    
    predictionTmp = wS.*deltaS + c;
    hold on;
    for j = 1:length(predictionTmp)
        if deltaRB1(j)<0 && deltaRB2(j)>0
            colorPlot = [0.56 0 0];
            plot(predictionTmp(j),meanChange(j),'o','LineWidth',1,'Color',colorPlot,'MarkerFaceColor',colorPlot);
        elseif deltaRB1(j)>0 && deltaRB2(j)<0
            colorPlot = [0.00 0 1];
            plot(predictionTmp(j),meanChange(j),'o','LineWidth',1,'Color',colorPlot,'MarkerFaceColor',colorPlot);
        else
            colorPlot = [0.56 0 1];
        end
        plot(predictionTmp(j),meanChange(j),'o','LineWidth',1,'Color',colorPlot,'MarkerFaceColor',colorPlot);
    end    
    plot([-3 3],[-3 3],'k--');
    xlim([-3 3]);
    ylim([-3 3]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaA');
    ylabel('Measured \DeltaA');
    title(['Correlation = ' num2str(corr(predictionTmp,meanChange),3)]);
    axis square;
end

