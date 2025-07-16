%% LOAD MAIN EXPERIMENT FILES

subjNum = 13;
bSave = false;
filePath = '/Users/benjaminchin/Documents/ARchromaScraps/meeting_Sept25/';

if subjNum==11
   % blockNums = 2:7;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   blockNums = 11:16;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};   
   subjName = ['S' num2str(subjNum) '-OD'];
elseif subjNum==12
   % blockNums = 2:7;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   % subjName = ['S' num2str(subjNum) '-OD'];

   blockNums = [11:15 18];
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];
elseif subjNum==13
   % blockNums = 3:8;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   blockNums = 12:17;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};   
   subjName = ['S' num2str(subjNum) '-OD'];   
elseif subjNum==14
   blockNums = 9:14;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   % blockNums = 3:8;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};   
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==15
   blockNums = 3:8;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];   
elseif subjNum==17
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==18
   blockNums = 10:15;
   % blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==19
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];         
elseif subjNum==20
   blockNums = 3:8;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];         
elseif subjNum==21
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];            
elseif subjNum==22
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];          
elseif subjNum==26
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==27
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==28
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];   
elseif subjNum==30
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
end

meanC = [];
rgb1all = [];
meanv00all = [];
nIndBadTracker = [];
c4all = {};
timeCell = {};

for i = 1:length(blockNums)
    blockNumTmp = blockNums(i);
    trialNumsTmp = trialNums{i};
    AFCp = ARCloadFileBVAMS(subjNum,blockNumTmp);
    for j = 1:length(trialNumsTmp)
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumsTmp(j),0);
        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65.
        indBadPupil = table2array(ZernikeTable(:,5))==0;
        PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(~indBadPupil,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0;
        nIndBadTracker(end+1) = sum(indBad);
        c(indBad,4) = mean(c(~indBad,4));
        if subjNum==18
           meanC(end+1,:) = mean(c(1:end,:),1);
        else
           meanC(end+1,:) = mean(c(1:end,:),1); % TAKE MEAN OF COEFFICIENTS    
        end
        c4all{end+1} = c(:,4);
        timeCell{end+1} = seconds(TimeStamp)-seconds(TimeStamp(1));
    end
    rgb1all = [rgb1all; AFCp.rgb100(trialNumsTmp,:)];
    meanv00all = [meanv00all; AFCp.meanv00(trialNumsTmp)./1.2255];
end

%% PLOTTING ALL TRIALS FOCUS FOR EACH OF 3 OPTICAL DISTANCES

defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
defocusAt550 = humanWaveDefocus(875)-humanWaveDefocus(550)+meanC(:,4)./defocusCorrectionFactor;
defocusAt875 = meanC(:,4)./defocusCorrectionFactor;

figure 
plot(meanv00all,defocusAt550,'ko');
set(gca,'FontSize',15);
xlim([1 4]);
hold on;
plot([1 4],[1 4],'k--','LineWidth',1);
axis square;
xlabel('Stimulus optical distance');
ylabel('Raw refractive power (D)');
title(['Subject ' num2str(subjNum-10)]);
if bSave
    saveas(gcf,[filePath 'stimRsp/S' num2str(subjNum) 'stimRsp'],'epsc');
end

%% PLOTTING ALL TRIAL MEANS PER CONDITION AND DISTANCE

lumScaleRGB = [4.0888 9.6669 1];

gammaRGB = [2.5 2.7 2.3];

rgbLumNorm = [lumScaleRGB(1).*rgb1all(:,1).^gammaRGB(1) lumScaleRGB(2).*rgb1all(:,2).^gammaRGB(2) lumScaleRGB(3).*rgb1all(:,3).^gammaRGB(3)];

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

for j = 1:length(optDistToCheckAll)
    figure;
    set(gcf,'Position',figPositions(j,:));
    hold on;
    optDistToCheck = optDistToCheckAll(j);
    indDist = abs(meanv00all-optDistToCheck)<0.01;
    for i = 1:size(conditionsOrderedNorm,1)
        ind = abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
              abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
              abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01 & ...
              abs(meanv00all-optDistToCheck)<0.01;
        if i<size(conditionsOrderedNorm,1)
            plot(i.*ones([sum(ind) 1]),defocusAt550(ind),'o','Color',conditionsOrderedNorm(i,:),'MarkerFaceColor',conditionsOrderedNorm(i,:));
        else
            plot(i.*ones([sum(ind) 1]),defocusAt550(ind),'o','Color','k','MarkerFaceColor','k');
        end
        defocusAt550mean(i) = mean(defocusAt550(ind));
    end
    plot(defocusAt550mean(1:5),'k-');
    plot(6:10,defocusAt550mean(6:10),'k-');
    plot([0 11],defocusAt550mean(11).*[1 1],'k--','LineWidth',1);
    xlim([0 11]);
    ylim(mean(defocusAt550(indDist))+[-0.6 0.6]);
    title(['Subject ' num2str(subjNum-10) ', Optical Distances = ' num2str(optDistToCheck)]);
    plot(5.5.*[1 1],ylim,'k-');
    set(gca,'FontSize',15);
    xlabel('Condition');
    ylabel('Defocus at 550nm');
    if bSave
       saveas(gcf,[filePath 'colorConditions/S' num2str(subjNum) 'ColorSplitDistInd' num2str(j)],'epsc');
    end
end

%% PLOTTING WAVELENGTH IN FOCUS

lumScaleRGB = [4.0888 9.6669 1];

gammaRGB = [2.5 2.7 2.3];

rgbLumNorm = [lumScaleRGB(1).*rgb1all(:,1).^gammaRGB(1) lumScaleRGB(2).*rgb1all(:,2).^gammaRGB(2) lumScaleRGB(3).*rgb1all(:,3).^gammaRGB(3)];

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

for j = 1:length(optDistToCheckAll)
    figure;
    set(gcf,'Position',figPositions(j,:));
    hold on;
    optDistToCheck = optDistToCheckAll(j);
    indDist = abs(meanv00all-optDistToCheck)<0.01;
    for i = 1:size(conditionsOrderedNorm,1)
        ind = abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
              abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
              abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01 & ...
              abs(meanv00all-optDistToCheck)<0.01;
        diffFromOptDist = defocusAt550(ind)-meanv00all(ind);
        wvInFocus = humanWaveDefocusInvert550anchor(diffFromOptDist);
        if i<size(conditionsOrderedNorm,1)
            plot(i.*ones([sum(ind) 1]),wvInFocus,'o','Color',conditionsOrderedNorm(i,:),'MarkerFaceColor',conditionsOrderedNorm(i,:));
        else
            plot(i.*ones([sum(ind) 1]),wvInFocus,'o','Color','k','MarkerFaceColor','k');
        end
        diffFromOptDistMean(i) = mean(diffFromOptDist);
        wvInFocusMean(i) = humanWaveDefocusInvert550anchor(diffFromOptDistMean(i));
    end
    plot(wvInFocusMean(1:5),'k-');
    plot(6:10,wvInFocusMean(6:10),'k-');
    plot([0 11],wvInFocusMean(11).*[1 1],'k--','LineWidth',1);
    xlim([0 11]);
    ylim([400 700]);
    title(['Subject ' num2str(subjNum-10) ', Optical Distances = ' num2str(optDistToCheck)]);
    plot(5.5.*[1 1],ylim,'k-');
    set(gca,'FontSize',15);
    xlabel('Condition');
    ylabel('Wavelength in focus (nm)');
    if bSave
       saveas(gcf,[filePath 'colorConditions/S' num2str(subjNum) 'ColorSplitDistInd' num2str(j)],'epsc');
    end
end


%% PLOT INDIVIDUAL TRIALS FOR REFERENCE

optDistToCheckAll = [1.5 2.5 3.5];
lengthTrialMax = 150;

for k = 1:length(optDistToCheckAll)
    optDistToCheck = optDistToCheckAll(k);
    indDist = abs(meanv00all-optDistToCheck)<0.01;
    figure;
    set(gcf,'Position',[148 265 1384 710]);
    hold on;
    for i = 1:(size(conditionsOrderedNorm,1)-1)
        subplot(2,5,i);
        hold on;
        ind = find(abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
                   abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
                   abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01 & ...
                   abs(meanv00all-optDistToCheck)<0.01);
        trialTmp875All = [];
        for j = 1:length(ind)
        % for j = 2
            trialTmp = zeros([1 lengthTrialMax]);
            trialTmp(1:length(c4all{ind(j)})) = c4all{ind(j)};
            trialTmp(trialTmp==0) = NaN;
            trialTmp550 = humanWaveDefocus(875)-humanWaveDefocus(550)+trialTmp./defocusCorrectionFactor;
            trialTmp875 = trialTmp./defocusCorrectionFactor;
            timeStampTmp = timeCell{ind(j)};
            plot(timeStampTmp,-trialTmp875(1:length(timeStampTmp)),'-','Color',conditionsOrderedNorm(i,:));
            trialTmp875All = [trialTmp875All; trialTmp875'];
        end
        % plot([timeStampTmp(1) timeStampTmp(end)],-mean(trialTmp875All(~isnan(trialTmp875All))).*[1 1],'-','Color',conditionsOrderedNorm(i,:),'LineWidth',2);
        errorbar(1.5,-mean(trialTmp875All(~isnan(trialTmp875All))),std(trialTmp875All(~isnan(trialTmp875All))), ...
                 'o','Color',conditionsOrderedNorm(i,:),'LineWidth',2, ...
                 'MarkerFaceColor','w','MarkerSize',10);
        axis square;
        xlim([0 3.1]);
        if subjNum==18
            ylim(mean(defocusAt875(indDist))+[-1.2 1.2]);
            yBuffer = 1.2;
        elseif subjNum==19
            ylim(mean(defocusAt875(indDist))+[-3 3]);  
            yBuffer = 3;
        elseif subjNum==17
            ylim(mean(defocusAt875(indDist))+[-1.2 1.2]);   
            yBuffer = 1.2;
        else
            ylim(mean(-defocusAt875(indDist))+[-0.6 0.6]);
            yBuffer = 0.6;
        end
        set(gca,'FontSize',15);
        trialStr = [];
        for l = 1:length(ind)
            trialStr = [trialStr ' ' num2str(ind(l))];
            text(50,mean(defocusAt875(indDist))+0.5,trialStr);
        end
        if i==1
            xlabel('Frame');
            ylabel('Defocus at 550nm');
            title(['Subject ' num2str(subjNum-10) ', Optical Distance = ' num2str(optDistToCheck)]);
        end
    end
    if bSave
       saveas(gcf,[filePath 'rawTraces/S' num2str(subjNum) 'TrialSplitDistInd' num2str(k)],'epsc');
    end
end

%% RSP AS FUNCTION OF LUMINANCE

subjNum = 12;

if subjNum==12
   blockNums = 8;
   trialNums = {1:30};
   subjName = ['S' num2str(subjNum) '-OD'];
elseif subjNum==13
   blockNums = 9;
   trialNums = {1:30};
   subjName = ['S' num2str(subjNum) '-OD'];   
end

meanC = [];
rgb1all = [];
meanv00all = [];
nIndBadTracker = [];
c4all = {};

for i = 1:length(blockNums)
    blockNumTmp = blockNums(i);
    trialNumsTmp = trialNums{i};
    AFCp = ARCloadFileBVAMS(subjNum,blockNumTmp);
    for j = 1:length(trialNumsTmp)
        [ZernikeTable, ~, ~, ~] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumsTmp(j),0);
        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65.
        indBadPupil = table2array(ZernikeTable(:,5))==0;
        PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(~indBadPupil,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0;
        nIndBadTracker(end+1) = sum(indBad);
        c(indBad,4) = mean(c(~indBad,4));
        meanC(end+1,:) = mean(c,1); % TAKE MEAN OF COEFFICIENTS    
        c4all{end+1} = c(:,4);
    end
    rgb1all = [rgb1all; AFCp.rgb100];
    meanv00all = [meanv00all; AFCp.meanv00./1.2255];
end

%% PLOT MEANS

defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
defocusAt550 = humanWaveDefocus(875)-humanWaveDefocus(550)+meanC(:,4)./defocusCorrectionFactor;
indHigherLum = abs(rgb1all(:,3)-1)<0.001;

figure; 
hold on;
plot(meanv00all(indHigherLum)-0.1,defocusAt550(indHigherLum),'ko','MarkerFaceColor','w');
plot(meanv00all(~indHigherLum)+0.1,defocusAt550(~indHigherLum),'ko','MarkerFaceColor','k');
set(gca,'FontSize',15);
xlim([1 4]);
axis square;
xlabel('Stimulus optical distance');
ylabel('Raw refractive power (D)');

%% LOAD RECORDINGS DURING ACUITY TASK

[ZernikeTable, ~, ~, ~] = ARCloadFileFIATallInstances('S12-OD',10,0);
NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65.
indBadPupil = table2array(ZernikeTable(:,5))==0;
PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
PARAMS.PupilFitSize=mean(table2array(ZernikeTable(~indBadPupil,5))); 
PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));

indBad = c(:,4)==0;
cGood = c(~indBad,4);
defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
defocusAt550 = humanWaveDefocus(875)-humanWaveDefocus(550)+cGood./defocusCorrectionFactor;

