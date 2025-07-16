%% ACTUAL POWER ANALYSIS

weightsR = [0.7125 0.1385 0.2997 0.2011 0.2584 0.2179 0.2376 0.2932];
weightsB = [-0.5044 -0.2884 -0.6663 -0.1809 -0.3804 -0.1181 -0.4334 -0.5439];
addVarExplClr = weightsR-weightsB;

%% JUST CHECKING EXPRESSION FOR STD ERROR OF THE MEAN

for i = 1:10000
    smpAddVar = randsample(addVarExplClr,8,true);
    muSmpAddVar(i) = mean(smpAddVar);
end

xSupport = -1:0.001:1;
pMu = normpdf(xSupport,0,std(muSmpAddVar));

figure;
plot(xSupport,pMu,'k','LineWidth',1.5); hold on;
plot(mean(muSmpAddVar).*[1 1],ylim,'k--','LineWidth',1.5);
set(gca,'FontSize',15);
xlabel('Sample mean'); ylabel('Probability');

%% COMPARE TWO-PARAMETER VS ONE-PARAMETER COLOR MODELS

twoParamR2 = [0.946 0.863 0.939 0.943 0.528];

oneParamR2 = [0.941 0.854 0.938 0.943 0.505];

figure;
hold on;
bar([1 3 5 7 9],twoParamR2,0.4,'FaceColor','w');
bar([2 4 6 8 10],oneParamR2,0.4,'FaceColor','k');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Correlation');
set(gca,'XTick',[1.5 3.5 5.5 7.5 9.5]);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5'});
legend({'2-parameter' '1-parameter'});

%%

varColor = [0.1989 0.0924 0.1340 0.1253 0.0462];

varAccStim = [0.6708 0.6084 0.7140 0.7344 0.2107];

figure;
hold on;
bar([1 3 5 7 9],varAccStim,0.4,'FaceColor','w');
bar([2 4 6 8 10],varColor,0.4,'FaceColor','m');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Variance explained');
set(gca,'XTick',[1.5 3.5 5.5 7.5 9.5]);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5'});
legend({'Accommodative demand' 'Color'});

%%

varColor = [0.1989 0.0924 0.1340 0.1253 0.0462];

varAccStim = [0.6708 0.6084 0.7140 0.7344 0.2107];

figure;
hold on;
bar([1 3 5 7 9],varAccStim,0.4,'FaceColor','w');
bar([2 4 6 8 10],addVarExplClr,0.4,'FaceColor','m');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Variance explained');
set(gca,'XTick',[1.5 3.5 5.5 7.5 9.5]);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5'});
legend({'Accommodative demand' 'Color'});

%% COMPARE WEIGHTING MODEL AND SWITCHING MODEL

% weightingModelRho = [0.946 0.863 0.939 0.943 0.528 0.971 0.853];
% switchingModelRho = [0.901 0.799 0.905 0.925 0.51  0.943 0.835];
% 
% 
% figure;
% hold on;
% bar([1 3 5 7 9 11 13],weightingModelRho.^2,0.4,'FaceColor','w');
% bar([2 4 6 8 10 12 14],switchingModelRho.^2,0.4,'FaceColor','k');
% set(gca,'FontSize',15);
% xlabel('Subject'); ylabel('Variance Explained');
% set(gca,'XTick',[1.5 3.5 5.5 7.5 9.5 11.5 13.5]);
% set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6' '7'});
% legend({'Weighting' 'Switching'});

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
    [weightsRBS1_x, weightsRBS1_y, rhoFull, rhoNoColor, rhoColor, aicLin, aicNoColor, weightsRBSci] = ARCnlz_linearModelnobias(sn(i),false,1000,0);
    rhoFullAll(i) = rhoFull;
    rhoNoColorAll(i) = rhoNoColor;
    weightsRBSall(i,:) = weightsRBS1_x;
    rhoColorAll(i) = rhoColor;
    aicLinAll(i) = aicLin;
    aicNoColorAll(i) = aicNoColor;
    weightsRBSciAll(:,:,i) = weightsRBSci;
    display(['Done with subject ' num2str(i)]);
end

%% PLOT WEIGHTS

labelsCell = {};

figure;
% set(gcf,'Position',[262 314 1239 594]);
for i = 1:size(weightsRBSall,1)
    weightsRBSciTmp = squeeze(weightsRBSciAll(:,:,i));
    hold on;
    bar((3*(i-1)+1),weightsRBSall(i,1),'FaceColor','r');
    plot((3*(i-1)+1).*[1 1],[weightsRBSciTmp(1,1) weightsRBSciTmp(2,1)],'k-');
    bar((3*(i-1)+2),weightsRBSall(i,2),'FaceColor','b');
    plot((3*(i-1)+2).*[1 1],[weightsRBSciTmp(1,2) weightsRBSciTmp(2,2)],'k-');
%     bar(3,weightsRBSall(i,3),'FaceColor','k');
%     plot([3 3],[weightsRBSciTmp(1,3) weightsRBSciTmp(2,3)],'k-');
    labelsCell{i} = ['S' num2str(sn(i))];
    set(gca,'FontSize',20);
    ylim(1.*[-2 2]);
    xlabel('Subject #');
    ylabel('Weight');
end
set(gca,'XTickLabel',labelsCell);
set(gca,'XTick',1.5:3:28.5);

figure;
% set(gcf,'Position',[262 314 1239 594]);
for i = 1:size(weightsRBSall,1)
    weightsRBSciTmp = squeeze(weightsRBSciAll(:,:,i));
    hold on;
    bar(i,weightsRBSall(i,3),'FaceColor','w','LineWidth',1);
    plot(i.*[1 1],[weightsRBSciTmp(1,3) weightsRBSciTmp(2,3)],'k-');
    set(gca,'XTick',[1:10]);
    set(gca,'XTickLabel',labelsCell);
    set(gca,'FontSize',20);
    ylim(1.*[0 1.15]);
    xlabel('Subject #');
    ylabel('Weight');
end

% Rscale = [1 1 1 1 1];
% figure;
% set(gca,'FontSize',15);
% hold on;
% plot([0 0.3],[0 1],'k--');
% for i = 1:length(dAll)
%     plot(dAll(i),Rscale(i)*weightsRBSall(i,1)-weightsRBSall(i,2),'ko','MarkerSize',10,'MarkerFaceColor','w');
%     text(dAll(i)+0.05,Rscale(i)*weightsRBSall(i,1)-weightsRBSall(i,2),['S' num2str(i)]);
% end
% axis square;
% xlabel('d'); ylabel('w_R - w_B');
% xlim([0 1.1]); ylim([0 1.1]);

%% PLOT VARIANCE EXPLAINED COMPARISONS

% figure;
% hold on;
% bar([1 4 7 10 13 16 19 22],rhoFullAll.^2,0.3,'FaceColor','w');
% % bar([2 5 8 11 14],rhoSwitchAll.^2,0.3,'FaceColor',[0.5 0.5 0.5]);
% bar([2 5 8 11 14 17 20 23],rhoNoColorAll.^2,0.3,'FaceColor','k');
% set(gca,'FontSize',15);
% xlabel('Subject'); ylabel('Proportion Variance Explained');
% set(gca,'XTick',[1.5 4.5 7.5 10.5 13.5 16.5 19.5 22.5]);
% set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6' '7' '8'});
% legend({'Weighting' 'No color'});

figure;
hold on;
bar(1:3:28,aicLinAll,0.3,'FaceColor','w');
% bar([2 5 8 11 14],aicSwitchAll,0.3,'FaceColor',[0.5 0.5 0.5]);
bar(2:3:29,aicNoColorAll,0.3,'FaceColor','k');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('AIC');
set(gca,'XTick',1.5:3:28.5);
set(gca,'XTickLabel',labelsCell);
legend({'Weighting' 'No color'});

%% PLOT VARIANCE EXPLAINED BY COLOR COMPARISONS

figure;
hold on;
bar([1 3 5 7 9 11],rhoColorAll.^2,0.4,'FaceColor','m');
bar([2 4 6 8 10 12],rhoColorSwitchAll.^2,0.4,'FaceColor',[0.5 0 0.5]);
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Variance Explained');
set(gca,'XTick',[1.5 3.5 5.5 7.5 9.5 11.5]);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6'});
legend({'Weighting' 'Switching'});

%% COMPARE WEIGHTS

redWeights = [0.7076 0.1043 0.3660 0.4499 0.2045 0.4052 0.2148];
blueWeights = [-0.2690 -0.4902 -0.2788 -0.3552 -0.0726 -0.2039 -0.1898];

%% COMPARE WITH AND WITHOUT BIAS TERM

sn = [17 11 16 19 21 23 24 25];
rhoFullAllnobias = [];
rhoFullBias = [];

for i = 1:length(sn)
    [weightsRBS1_xNoBias, weightsRBS1_yNoBias, rhoFull, rhoNoColor, rhoColor] = ARCnlz_linearModelnobias(sn(i),false);
    rhoFullAllnobias(i) = rhoFull;
    display(['Done with subject ' num2str(i)]);
end

for i = 1:length(sn)
    [weightsRBS1_x, weightsRBS1_y, rhoFull, rhoNoColor, rhoColor] = ARCnlz_linearModel(sn(i),false);
    rhoFullAllBias(i) = rhoFull;
    display(['Done with subject ' num2str(i)]);
end

%% PLOT VARIANCE EXPLAINED WITH AND WITHOUT BIAS TERM

figure;
hold on;
bar([1 3 5 7 9 11 13 15],rhoFullAllBias.^2,0.4,'FaceColor','k');
bar([2 4 6 8 10 12 14 16],rhoFullAllnobias.^2,0.4,'FaceColor','w');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Variance Explained');
set(gca,'XTick',[1.5 3.5 5.5 7.5 9.5 11.5 13.5 15.5]);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6' '7' '8'});
legend({'With bias' 'No bias'});

%% COMPARE 3-PARAMETER AND 2-PARAMETER MODELS

sn = [18 23 24 26 27];
rhoFullAll = [];
rhoNoColorAll = [];
weightsRBSall = [];
rhoColorAll = [];

for i = 1:length(sn)
    [weightsRBS1_x, weightsRBS1_y, rhoFull, rhoNoColor, rhoColor] = ARCnlz_linearModelnobias(sn(i),false,1000);
    rhoFullAll(i) = rhoFull;
    rhoNoColorAll(i) = rhoNoColor;
    weightsRBSall(i,:) = weightsRBS1_x;
    rhoColorAll(i) = rhoColor;
    display(['Done with subject ' num2str(i)]);
end

save('/home/ben/Documents/ARchroma/params2paramModel','rhoFullAll','rhoNoColorAll','weightsRBSall');

%% COMPARE PARAMETERS FROM ALL TRIALS VS RANDOM ONLY TRIALS

load('/Users/benchin/Documents/ARchroma/cache/paramsAll.mat');
rhoFullAllAll = rhoFullAll;
weightsRBSallAll = weightsRBSall;
load('/Users/benchin/Documents/ARchroma/cache/paramsRandom.mat');
rhoFullAllRandom = rhoFullAll;
weightsRBSallRandom = weightsRBSall;

figure;
set(gcf,'Position',[318 247 1135 708]);
subplot(2,2,1);
set(gca,'FontSize',15);
hold on;
plot(rhoFullAllAll,rhoFullAllRandom,'ko','MarkerSize',10);
plot([0 1],[0 1],'k--');
xlabel('All trials'); ylabel('Random trials only'); title('Correlation');
axis square;
subplot(2,2,2);
set(gca,'FontSize',15);
hold on;
plot(weightsRBSallAll(:,3),weightsRBSallRandom(:,3),'ko','MarkerSize',10);
plot([0 1],[0 1],'k--');
xlabel('All trials'); ylabel('Random trials only'); title('Gain');
axis square;
subplot(2,2,3);
set(gca,'FontSize',15);
hold on;
plot([0 1],[0 1],'k--');
plot(weightsRBSallAll(:,1),weightsRBSallRandom(:,1),'ro','MarkerSize',10,'MarkerFaceColor','w');
xlabel('All trials'); ylabel('Random trials only'); title('Red weight');
axis square;
subplot(2,2,4);
set(gca,'FontSize',15);
hold on;
plot([-1 0],[-1 0],'k--');
plot(weightsRBSallAll(:,2),weightsRBSallRandom(:,2),'bo','MarkerSize',10,'MarkerFaceColor','w');
xlabel('All trials'); ylabel('Random trials only'); title('Blue weight');
axis square;

%% COMPARE 3-PARAMETER MODEL WITH 2-PARAMETER MODEL

load('/home/ben/Documents/ARchroma/params3paramModel.mat');
rhoFullAll3param = rhoFullAll;
weightsRBSall3param = weightsRBSall;
rhoNoColorAll3param = rhoNoColorAll;
load('/home/ben/Documents/ARchroma/params2paramModel.mat');
rhoFullAll2param = rhoFullAll;
weightsRBSall2param = weightsRBSall;
rhoNoColorAll2param = rhoNoColorAll;

figure;
hold on;
bar([1 4 7 10 13],rhoFullAll3param.^2,0.3,'FaceColor','w');
bar([2 5 8 11 14],rhoFullAll2param.^2,0.3,'FaceColor',[0.7 0.7 0.7]);
bar([3 6 9 12 15],rhoNoColorAll2param.^2,0.3,'FaceColor',[0 0 0]);
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Variance explained');
set(gca,'XTick',[2 5 8 11 14]);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5'});
legend({'3-parameter' '2-parameter' 'No color'});

figure;
set(gcf,'Position',[251 387 1120 420]);
set(gca,'FontSize',15);
subplot(1,2,1);
hold on;
plot(weightsRBSall2param(:,1),weightsRBSall3param(:,1),'ro','MarkerSize',10);
plot([0 1],[0 1],'k--');
xlabel('2-param ratio'); ylabel('3-param w_R');
axis square;
subplot(1,2,2);
hold on;
plot(weightsRBSall2param(:,1),weightsRBSall3param(:,2),'bo','MarkerSize',10);
plot([-0.1 1],[0.1 -1],'k--');
xlabel('2-param ratio'); ylabel('3-param w_B');
axis square;

figure;
hold on;
bar([1 3 5 7 9 11 13 15],rhoFullAll2param.^2,0.4,'FaceColor','w');
bar([2 4 6 8 10 12 14 16],rhoNoColorAll2param.^2,0.4,'FaceColor','k');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Variance explained');
set(gca,'XTick',[1.5 3.5 5.5 7.5 9.5 11.5 13.5 15.5]);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6' '7' '8'});
legend({'With color' 'No color'});

figure;
bar(1:8,rhoFullAll2param.^2 - rhoNoColorAll2param.^2,0.4,'FaceColor','w');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Variance explained due to color');
set(gca,'XTick',1:9);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6' '7' '8'});

%% AIC ANALYSIS

sn = [18 23 24 26 27 28 29 30 32 33];
dAll1 = [];
rbThreshold1 = [];
aicSwitchAll1 = [];
dAll2 = [];
rbThreshold2 = [];
aicSwitchAll2 = [];
aicLinAll = [];
aicNoColorAll = [];
% meanDiffPureRedVsBlue = [0.9607 1.0299 0.9122 0.1594 0.4769 0.5201 0.2107 0.2619 0.8158 0.3803];
% meanDiffPureRedVsBlue = [0.8154 0.7792 0.5180 0.6472 0.4920 0.3267 0.5820 0.5690 0.7679 0.6513];
meanDiffPureRedVsBlue = 0.8.*ones([1 10]);

for i = 1:length(sn)
    [d, wS, rbThreshold, rhoSwitch, rhoColorSwitch, aicSwitch] = ARCnlz_linearSwitching(sn(i),false,0);
    dAll1(i) = d;
    rbThreshold1(i) = rbThreshold;
    aicSwitchAll1(i) = aicSwitch;
    [d, wS, rbThreshold, rhoSwitch, rhoColorSwitch, aicSwitch] = ARCnlz_linearSwitching(sn(i),false,meanDiffPureRedVsBlue(i));
    dAll2(i) = d;
    rbThreshold2(i) = rbThreshold;
    aicSwitchAll2(i) = aicSwitch;    
    display(['Done with subject ' num2str(i)]);
end

for i = 1:length(sn)
    [weightsRBS1_x, weightsRBS1_y, rhoFull, rhoNoColor, rhoColor, aicLin, aicNoColor, weightsRBSci] = ARCnlz_linearModelnobias(sn(i),false,1000,0);
    aicLinAll(i) = aicLin;
    aicNoColorAll(i) = aicNoColor;
    display(['Done with subject ' num2str(i)]);
end

%% PLOT AIC ANALYSIS

figure;
set(gcf,'Position',[215 337 792 420]);
hold on;
bar([1 6 11 16 21 26 31 36 41 46],aicLinAll,0.18,'FaceColor','w');
bar([2 7 12 17 22 27 32 37 42 47],aicSwitchAll1,0.18,'FaceColor',[0.5 0.5 0.5]);
bar([3 8 13 18 23 28 33 38 43 48],aicSwitchAll2,0.18,'FaceColor',[0.25 0.25 0.25]);
bar([4 9 14 19 24 29 34 39 44 49],aicNoColorAll,0.18,'FaceColor','k');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('AIC');
set(gca,'XTick',[2.5 7.5 12.5 17.5 22.5 27.5 32.5 37.5 42.5 47.5]);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'});
legend({'Weighting' 'Switch Float d' 'Switch fixed d' 'No color'});

figure;
set(gcf,'Position',[215 337 692 420]);
hold on;
bar([1 5 9  13 17 21 25 29 33 37],aicLinAll,0.18,'FaceColor','w');
bar([2 6 10 14 18 22 26 30 34 38],aicSwitchAll2,0.18,'FaceColor',[0.5 0.5 0.5]);
bar([3 7 11 15 19 23 27 31 35 39],aicNoColorAll,0.18,'FaceColor','k');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('AIC');
set(gca,'XTick',[2:4:38]);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'});
legend({'Weighting' 'Switching' 'No color'});

figure;
set(gcf,'Position',[215 337 792 420]);
hold on;
plot([0 170],[0 170],'k--');
plot(aicNoColorAll,aicLinAll,'ko','MarkerFaceColor',1.*[1 1 1],'MarkerSize',15);
xlim([0 170]);
ylim([0 170]);
set(gca,'FontSize',15);
axis square;
xlabel('AIC (no color model)'); ylabel('AIC (linear weighting model)');

figure;
set(gcf,'Position',[215 337 792 420]);
hold on;
plot([0 170],[0 170],'k--');
plot(aicSwitchAll2,aicLinAll,'ko','MarkerFaceColor',0.5.*[1 1 1],'MarkerSize',15);
xlim([0 170]);
ylim([0 170]);
set(gca,'FontSize',15);
axis square;
xlabel('AIC (switching model)'); ylabel('AIC (linear weighting model)');

figure;
set(gcf,'Position',[215 337 792 420]);
hold on;
plot([0 170],[0 170],'k--');
plot(aicNoColorAll,aicSwitchAll2,'ko','MarkerFaceColor',0.5.*[1 1 1],'MarkerSize',15);
xlim([0 170]);
ylim([0 170]);
set(gca,'FontSize',15);
axis square;
xlabel('AIC (no color model)'); ylabel('AIC (switching model)');

figure;
set(gcf,'Position',[215 337 792 420]);
hold on;
% bar([1 3 5 7 9  11 13 15 17 19],dAll1,0.4,'FaceColor','w');
% bar([2 4 6 8 10 12 14 16 18 20],dAll2,0.4,'FaceColor',[0.5 0.5 0.5]);
bar(1:10,dAll2,0.4,'FaceColor',[0.5 0.5 0.5]);
% plot(xlim,0.35.*[1 1],'k--');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Switching d');
% set(gca,'XTick',[1.5 3.5 5.5 7.5 9.5 11.5 13.5 15.5 17.5 19.5]);
set(gca,'XTick',1:10);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'});
% legend({'Switch Float d' 'Switch constrained d'});

figure;
set(gcf,'Position',[215 337 792 420]);
hold on;
% bar([1 3 5 7 9  11 13 15 17 19],rbThreshold2,0.4,'FaceColor','w');
% bar([2 4 6 8 10 12 14 16 18 20],rbThreshold2,0.4,'FaceColor',[0.5 0.5 0.5]);
bar(1:10,rbThreshold1,0.4,'FaceColor','w');
set(gca,'FontSize',15);
xlabel('Subject'); ylabel('Switching threshold');
set(gca,'XTick',1:10);
set(gca,'XTickLabel',{'1' '2' '3' '4' '5' '6' '7' '8' '9' '10'});
% legend({'Switch Float d' 'Switch 0.35 d'});
ylim([-0.5 0.5]);

%% COMPARING OUTPUTS OF LINEAR AND SWITCHING MODELS

subjNumber = [18 23 24 26 27 29 32 33];

figure;
set(gcf,'Position',[105 245 1479 647]);  
for i = 1:length(subjNumber)
    [d, wS, rbThreshold, rhoSwitch, rhoColorSwitch, aic, trialMeansSwitch] = ARCnlz_linearSwitching(subjNumber(i),0,0);
    [weightsRBS1_x, weightsRBS1_y, rhoFull, rhoNoColor, rhoColor, aic, aicNoColor, weightsRBSci, trialMeansLin] = ARCnlz_linearModelnobias(subjNumber(i),0,100,0);
    subplot(2,4,i);
    plot(trialMeansLin,trialMeansSwitch,'ko'); hold on;
    plotExtent = max([xlim ylim]);
    plot(plotExtent.*[-1 1],plotExtent.*[-1 1],'k--');
    xlim(plotExtent.*[-1 1]);
    ylim(plotExtent.*[-1 1]);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Linear model prediction');
    ylabel('Switching model prediction');    
end

