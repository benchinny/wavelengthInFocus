%% LOADING DATA

subjNum = 10;
objFunc = 'NLL';

if subjNum==10
    subjName = 'S20-OD';
    blockNumAll = 3:8;
elseif subjNum==3
    subjName = 'S13-OD';
    blockNumAll = 12:17;
elseif subjNum==1
    subjName = 'S11-OD';
    blockNumAll = 11:16;
elseif subjNum==5
    subjName = 'S15-OD';
    blockNumAll = 3:8;
elseif subjNum==9
    subjName = 'S19-OD';
    blockNumAll = 2:7;
elseif subjNum==16
    subjName = 'S26-OD';
    blockNumAll = 2:7;
elseif subjNum==17
    subjName = 'S27-OD';
    blockNumAll = 2:7;
elseif subjNum==18
    subjName = 'S28-OD';
    blockNumAll = 2:7;
elseif subjNum==20
    subjName = 'S30-OD';
    blockNumAll = 2:7;
end

trialNumAll = 1:36;

defocus875 = [];
optDistAll = [];
rgbAll = [];

for k = 1:length(blockNumAll)
    AFCp = ARCloadFileBVAMS(subjNum+10,blockNumAll(k));
    optDistAll = [optDistAll; AFCp.meanv00./1.2255];
    rgbAll = [rgbAll; AFCp.rgb100];
    for l = 1:length(trialNumAll)
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumAll(k),trialNumAll(l),0);

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS = struct;
        PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(:,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0;
        meanC = mean(c(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS  
        defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
        defocus875(end+1,:) = meanC(4)./defocusCorrectionFactor;
    end
end

rgbUnq = unique(rgbAll,'rows');

%% RMSE FOR DIFFERENT WEIGHTS

wLM = 0:0.1:1;
wS = 1-wLM;
RMSE = zeros(size(wLM));

parfor i = 1:length(wLM)
    RMSE(i) = ARCtestWvInFocusMeanZspatFilterObjFunc(subjNum,defocus875,rgbAll,optDistAll,[wLM(i).*[1 1] -wS(i)]);
    display(['Weight index' num2str(i)]);
end

%% PLOTTING WEIGHTS VS RMSE

figure;
hold on;
plot(wLM(~isnan(RMSE)),RMSE(~isnan(RMSE)),'k-','LineWidth',1);
plot(xlim,RMSEflat.*[1 1],'k--');
axis square;
set(gca,'FontSize',15);
ylabel('RMSE');
xlabel('Weight on L+M');

%% VISUALIZING PREDICTIONS VS ACTUAL DEFOCUS FOR DIFFERENT WEIGHTS

wLM = 0:0.1:1;
wS = 1-wLM;
coneWeightsFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsErrorSpatFilter/colorMechPredictions/';

rgbLumNorm = [];
rgbLumNorm(:,1) = (rgbUnq(:,1).^2.5)./0.2442;
rgbLumNorm(:,2) = (rgbUnq(:,2).^2.7)./0.1037;
rgbLumNorm(:,3) = (rgbUnq(:,3).^2.3)./1;
rgbLumNorm(rgbLumNorm>1) = 1;

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

for i = 1:size(conditionsOrderedNorm,1)
    ind(i) = find(abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
                  abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
                  abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01);
end

parfor k = 1:length(wLM)
    [~, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCtestWvInFocusMeanZspatFilterPlotHelper(subjNum,defocus875,rgbAll,optDistAll,[wLM(k).*[1 1] -wS(k)]);
    optDistTag = imresize(optDistUnq',size(defocus875mean),'nearest');
    [pFit,RMSE(k)] = ARCfitLagLead(defocus875predTmp(:),defocus875mean(:),optDistTag(:),true,objFunc);
    
    if k==1
       [pFitFlat,RMSEflat(k)] = ARCfitLagLead(optDistTag(:),defocus875mean(:),optDistTag(:),true,objFunc);
    end

    figure;
    set(gcf,'Position',[118 470 1453 420]);
    for i = 1:size(defocus875predTmp,2)
        subplot(1,3,i);
        hold on;
        plot(1:length(ind),defocus875predTmp(ind,i)-optDistUnq(i)*pFit(1)-pFit(2),'b-');
        plot(1:length(ind),defocus875mean(ind,i),'k-');
        for j = 1:length(ind)
            plot(j,defocus875mean(ind(j),i),'ko','MarkerFaceColor',conditionsOrderedNorm(j,:), ...
                 'MarkerSize',10);
        end
        axis square;
        set(gca,'FontSize',15);
        set(gca,'XTick',[]);
        xlabel('Condition');
        ylabel('Defocus at 875nm (D)');
        if i==1
            title(['Subject ' num2str(subjNum) ', Distance = ' num2str(optDistUnq(i))]);
        elseif i==2
            title(['Weights = [' num2str(wLM(k)) ' ' num2str(wLM(k)) ' ' num2str(wS(k)) '], Distance = ' num2str(optDistUnq(i))]);
        else
            title(['Distance = ' num2str(optDistUnq(i))]);
        end
    end
    saveas(gcf,[coneWeightsFolder 'LplusMminusSpredCont' num2str(subjNum) 'weight' num2str(k)],'png');
    display(['Cone weights ' num2str(k)]);
end

%% SEARCH INDIVIDUAL CONE WEIGHTS

wLM = 0.4:0.1:1.4;
wLprop = 0.3:0.05:0.75;
wS = -0.5;
coneWeightsFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsErrorSpatFilter/colorMechPredictions/';

rgbLumNorm = [];
rgbLumNorm(:,1) = (rgbUnq(:,1).^2.5)./0.2442;
rgbLumNorm(:,2) = (rgbUnq(:,2).^2.7)./0.1037;
rgbLumNorm(:,3) = (rgbUnq(:,3).^2.3)./1;
rgbLumNorm(rgbLumNorm>1) = 1;

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

for i = 1:size(conditionsOrderedNorm,1)
    ind(i) = find(abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
                  abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
                  abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01);
end

RMSEall = zeros([length(wLM) length(wLprop)]);

for l = 1:length(wLM)
    parfor k = 1:length(wLprop)
        wL = wLM(l)*wLprop(k);
        wM = wLM(l)-wL;
        [~, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCtestWvInFocusMeanZspatFilterPlotHelper(subjNum,defocus875,rgbAll,optDistAll,[wL wM wS]);
        optDistTag = imresize(optDistUnq',size(defocus875mean),'nearest');
        [pFit,RMSE(k)] = ARCfitLagLead(defocus875predTmp(:),defocus875mean(:),optDistTag(:),true,objFunc);
        if k==1
            [pFitFlat,RMSEflat(k)] = ARCfitLagLead(optDistTag(:),defocus875mean(:),optDistTag(:),true,objFunc);
        end
        
        display(['Weights = [' num2str(wL) ' ' num2str(wM) ' ' num2str(wS)]);
    end
    RMSEall(l,:) = RMSE;
end

save([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResults' num2str(round(-wS*10)) '.mat'],'RMSEall','wS','wLM','wLprop');

%% VISUALIZE BEST WEIGHTS

wL = 0.675;
wM = 0.225;
wS = -1;
coneWeightsFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsErrorSpatFilter/colorMechPredictions/';

rgbLumNorm = [];
rgbLumNorm(:,1) = (rgbUnq(:,1).^2.5)./0.2442;
rgbLumNorm(:,2) = (rgbUnq(:,2).^2.7)./0.1037;
rgbLumNorm(:,3) = (rgbUnq(:,3).^2.3)./1;
rgbLumNorm(rgbLumNorm>1) = 1;

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

for i = 1:size(conditionsOrderedNorm,1)
    ind(i) = find(abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
                  abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
                  abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01);
end

RMSEall = zeros([length(wL) length(wM)]);

for l = 1:length(wL)
    for k = 1:length(wM)
        [~, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCtestWvInFocusMeanZspatFilterPlotHelper(subjNum,defocus875,rgbAll,optDistAll,[wL(l) wM(k) wS]);
        optDistTag = imresize(optDistUnq',size(defocus875mean),'nearest');
        [pFit,RMSE(k)] = ARCfitLagLead(defocus875predTmp(:),defocus875mean(:),optDistTag(:),true,objFunc);
        
        if k==1
            [pFitFlat,RMSEflat(k)] = ARCfitLagLead(optDistTag(:),defocus875mean(:),optDistTag(:),true,objFunc);
        end

        figure;
        set(gcf,'Position',[118 470 1453 420]);
        for i = 1:size(defocus875predTmp,2)
            subplot(1,3,i);
            hold on;
            plot(1:length(ind),defocus875predTmp(ind,i)-optDistUnq(i)*pFit(1)-pFit(2),'k-');
            % plot(1:length(ind),defocus875mean(ind,i),'k-');
            for j = 1:length(ind)
                plot(j,defocus875mean(ind(j),i),'ko','MarkerFaceColor',conditionsOrderedNorm(j,:), ...
                     'MarkerSize',10);
            end
            axis square;
            set(gca,'FontSize',15);
            set(gca,'XTick',[]);
            xlabel('Condition');
            ylabel('Defocus at 875nm (D)');
            if i==1
                title(['Subject ' num2str(subjNum) ', Distance = ' num2str(optDistUnq(i))]);
            elseif i==2
                title(['Weights = [' num2str(wL(l)) ' ' num2str(wM(k)) ' ' num2str(wS) '], Distance = ' num2str(optDistUnq(i))]);
            else
                title(['Distance = ' num2str(optDistUnq(i))]);
            end
        end
        % saveas(gcf,[coneWeightsFolder 'LplusMminusSpredCont' num2str(subjNum) 'weight' num2str(k)],'png');
        display(['Weights = [' num2str(wL(l)) ' ' num2str(wM(k)) ' ' num2str(wS)]);
    end
end
