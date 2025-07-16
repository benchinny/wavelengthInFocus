%% LOAD DATA

subjNum = 20;

% PARAMETERS OF WAVEFRONT ANALYSIS
PARAMS.PupilSize = 7; %default values - will be replaced depending on choices below
PARAMS.PupilFieldSize =6; %default values - will be replaced depending on choices below
PARAMS.PupilFitSize = 7; %default values - will be replaced depending on choices below

wvfFiles = ARCacuAnalysisWvfSubj(subjNum);

dataFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/csvFiles/SUBJ/';

cAll = [];

for i = 1:length(wvfFiles)
    ZernikeTable = readtable([dataFolder wvfFiles{i}]);
    NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
    c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
    PARAMS = struct;
    PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
    PARAMS.PupilFitSize=mean(table2array(ZernikeTable(:,5))); 
    PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
    c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
    cAll = [cAll; c];
end

indBad = cAll(:,4)==0 | cAll(:,4)<-10 | cAll(:,4)>2.5;
meanC = mean(cAll(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS

defocusScaleFactor = 0.5774;

defocusOrig = meanC(4);
defocusOrigScaled = defocusOrig/defocusScaleFactor;

[unqFocDst,PC,PCci,dprime,dprimeCI,PCfit,dprimeFitAll,PCfitSupport] = ARCacuAnalysisSubjective(subjNum,0);
scaleFac = 0.8;

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/meeting_Sept25/allExp1DataRGB.mat');
indAcc = abs(optDistCndAll+2.5)<0.001 & abs(rgbLumNormCndAll(:,1)-1)<0.001 & ...
         abs(rgbLumNormCndAll(:,3)-1)<0.001 & abs(rgbLumNormCndAll(:,2)-0)<0.001 & ...
         subjNumTag==subjNum;
defocus875acc = mean(defocusAt875cellAll{indAcc});
diff875 = -defocus875acc-defocusOrigScaled;

%% PLOT ACUITY PERFORMANCE VS DISTANCE INCREMENT

figure;
hold on;
% plot(unqFocDst.*scaleFac,PC,'o-','Color',rgbUnq,'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
plot(PCfitSupport,PCfit,'-','Color',[0.56 0 1],'LineWidth',1.5);
plot(PCfitSupport+diff875,PCfit,'--','Color',[0.56 0 1],'LineWidth',1.5);
errorbar(unqFocDst.*scaleFac,PC,PC-PCci(1,:),PCci(2,:)-PC,'o','Color',[0.56 0 1],'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
axis square;
ylim([0.4 1]);
xlim([-1.2 1.2]);
formatFigure('Relative optical distance (D)','Proportion Correct',['Subject ' num2str(subjNum) ', D_{acu} = ' num2str(defocusOrigScaled,3) ', D_{acc} = ' num2str(-defocus875acc,3)]);
legend('Acuity','Accommodation');

[~,indPeak] = max(PCfit);
acuPeak = PCfitSupport(indPeak)+diff875;

%% ACUITY VS DISTANCE INCREMENT IN TERMS OF D-PRIME

figure;
hold on;
% plot(unqFocDst.*scaleFac,PC,'o-','Color',rgbUnq,'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
plot(PCfitSupport,dprimeFitAll,'-','Color',[0.56 0 1],'LineWidth',1.5);
errorbar(unqFocDst.*scaleFac,dprime,dprimeCI-dprimeCI(1,:),dprimeCI(2,:)-dprime,'o','Color',[0.56 0 1],'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
axis square;
% ylim([0.4 1]);
xlim([-1.2 1.2]);
formatFigure('Relative optical distance (D)','D-prime',['Subject ' num2str(subjNum) ', D_{acu} = ' num2str(defocusOrigScaled,3) ', D_{acc} = ' num2str(-defocus875acc,3)]);

dprimeMax = max(dprimeFitAll);
dprime0 = interp1(PCfitSupport,dprimeFitAll,0);
dprimeRatio = dprime0/dprimeMax;

%%

dprimeRatioAll = [2.41 1.000 1.000 1.41 1.000 1.98 1.76 1.023];

figure;
% boxplot(dprimeRatioAll,'LineWidth',1,'FaceColor',[0.56 0 1]);
boxplot(dprimeRatioAll);
set(gca,'FontSize',15);
hold on;
for i = 1:length(dprimeRatioAll)
    if abs(dprimeRatioAll(i)-1)<0.01
        plot(1+i*0.1,dprimeRatioAll(i),'ko','MarkerSize',10,'MarkerFaceColor',[0.56 0 1]);
    else
        plot(1,dprimeRatioAll(i),'ko','MarkerSize',10,'MarkerFaceColor',[0.56 0 1]);
    end
end
% set(gca,'XTick',1:8);
% set(gca,'XTickLabel',{'S1' 'S3' 'S5' 'S10' 'S16' 'S17' 'S18' 'S20'});
% xlabel('Subject');
ylabel('d''_{max}/d''_{0}');
% ylim([0 1]);

%% COMPARE ACCOMMODATION AND ACUITY EXPERIMENT FOR PURPLE 2.5D CONDITION

modelPredictionPurpleAt2pt5 = [1.36 1.756 1.864 1.633 1.463 1.815 1.355 1.603];
defocus875accAll = [1.138 1.659 1.828 1.522 1.294 1.629 1.155 1.604];
defocusOrigScaledAll = [1.224 1.4785 1.8230 1.4982 1.6005 1.5845 1.3686 1.6974];
defocusOrigScaledCIall = [1.1022 1.3805; 1.3391 1.6280; 1.7148 1.9500; 1.3001 1.6979; ...
                          1.4398 1.7393; 1.3628 1.7745; 1.1458 1.5710; 1.5912 1.8037];
defocusOrigScaledCI95 = [0.8100 1.5203; 1.1545 1.7632; 1.5140 2.0659; 1.1073 1.8363; ...
                         1.2657 1.8167; 1.1812 2.1370; 1.0553 1.7222; 1.4949 1.9042];

plotIndModelPredictionPurpleAt2pt5 = 1:4:29;
plotIndDefocus875accAll = 2:4:30;
plotIndDefocusOrigScaled = 3:4:31;

figure;
hold on;
plot(plotIndModelPredictionPurpleAt2pt5,modelPredictionPurpleAt2pt5,'o','MarkerFaceColor','w','MarkerSize',15,'LineWidth',1);
plot(plotIndDefocus875accAll,defocus875accAll,'o','MarkerFaceColor','w','MarkerSize',15,'LineWidth',1);
errorbar(plotIndDefocusOrigScaled,defocusOrigScaledAll,defocusOrigScaledAll-defocusOrigScaledCI95(:,1)',defocusOrigScaledCI95(:,2)'-defocusOrigScaledAll,'o','MarkerFaceColor','w','MarkerSize',15,'LineWidth',1);
set(gca,'FontSize',15);
ylim([0 2.5]);
ylabel('Defocus at 875nm');
xlabel('Subject');
set(gca,'XTick',2:4:30);
set(gca,'XTickLabel',{'S1' 'S3' 'S5' 'S10' 'S16' 'S17' 'S18' 'S20'});
legend('Model fit (acc exp)','Accommodation exp','Acuity exp');

figure;
hold on;
plot(plotIndModelPredictionPurpleAt2pt5,modelPredictionPurpleAt2pt5,'o','MarkerFaceColor','w','MarkerSize',15,'LineWidth',1);
plot(plotIndDefocus875accAll,defocus875accAll,'o','MarkerFaceColor','w','MarkerSize',15,'LineWidth',1);
errorbar(plotIndDefocusOrigScaled,defocusOrigScaledAll,defocusOrigScaledAll-defocusOrigScaledCIall(:,1)',defocusOrigScaledCIall(:,2)'-defocusOrigScaledAll,'o','MarkerFaceColor','w','MarkerSize',15,'LineWidth',1);
set(gca,'FontSize',15);
ylim([0 2.5]);
ylabel('Defocus at 875nm');
xlabel('Subject');
set(gca,'XTick',2:4:30);
set(gca,'XTickLabel',{'S1' 'S3' 'S5' 'S10' 'S16' 'S17' 'S18' 'S20'});
legend('Model fit (acc exp)','Accommodation exp','Acuity exp');

%% DISTANCE OF PEAK ACUITY PLOT

peakAcuityAll = [-0.7102 -0.4229 0.1432 -0.3498 -0.1546 -0.6124 -0.8071 -0.2844];

figure;
bar(peakAcuityAll,'LineWidth',1,'FaceColor',[0.56 0 1]);
set(gca,'FontSize',15);
set(gca,'XTick',1:8);
set(gca,'XTickLabel',{'S1' 'S3' 'S5' 'S10' 'S16' 'S17' 'S18' 'S20'});
xlabel('Subject');
ylabel('Best acuity relative distance (D)');
ylim([-1 1]);

