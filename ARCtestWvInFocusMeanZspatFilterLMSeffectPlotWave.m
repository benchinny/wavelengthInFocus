function [aic, pFit, wvMeanAll, wvPredAll] = ARCtestWvInFocusMeanZspatFilterLMSeffectPlotWave(subjNum,modelType)

coneWeightsFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsErrorSpatFilter/colorMechPredictions/';
objFunc = 'RMS';
modelCompFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/';

if strcmp(modelType,'LMS')
    wS = -1;
    if subjNum==20
        wS = -0.25;
    end
    if subjNum==5
        wS = -0.5;
    end
    load([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResults' num2str(round(-wS*10)) '.mat'],'RMSEall','wS','wLM','wLprop');
    
    [wLpropGrid,wLMgrid] = meshgrid(wLprop,wLM);
    
    indMin = RMSEall==min(RMSEall(:));
    wLMmin = wLMgrid(indMin);
    wLpropMin = wLpropGrid(indMin);
    wL = wLMmin*wLpropMin;
    wM = wLMmin-wL;
    nParams = 4;
end

if strcmp(modelType,'LM')
    wS = 0;
    load([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResults' num2str(round(-wS*10)) '.mat'],'RMSEall','wS','wLM','wLprop');
    
    [wLpropGrid,wLMgrid] = meshgrid(wLprop,wLM);
    
    indMin = RMSEall==min(RMSEall(:));
    wLMmin = wLMgrid(indMin);
    wLpropMin = wLpropGrid(indMin);
    wL = wLMmin*wLpropMin;
    wM = wLMmin-wL;
    nParams = 3;
end

if strcmp(modelType,'LminusM')
    wS = 0;
    load([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResultsLminusM.mat'],'RMSEall','wLM','wLprop');
    
    [wLpropGrid,wLMgrid] = meshgrid(wLprop,wLM);
    
    indMin = RMSEall==min(RMSEall(:));
    wLMmin = wLMgrid(indMin);
    wLpropMin = wLpropGrid(indMin);
    wL = wLMmin*wLpropMin;
    wM = -(wLMmin-wL);    
    nParams = 3;
end

if strcmp(modelType,'Lum')
    wL = 0.72;
    wM = 0.28;
    wS = 0;
    nParams = 2;
end

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

diffFromOptDist = defocus875-optDistAll;
% EXCLUDE DATA FOR WHICH PARTICIPANT WAS ACCOMMODATING OUTSIDE OF
% VISIBLE RANGE
indGood = humanWaveDefocusInvertARC(875,diffFromOptDist,subjNum)>380 & ...
          humanWaveDefocusInvertARC(875,diffFromOptDist,subjNum)<780;
defocus875 = defocus875(indGood);
rgbAll = rgbAll(indGood,:);
optDistAll = optDistAll(indGood);

for i = 1:size(conditionsOrderedNorm,1)
    ind(i) = find(abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
                  abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
                  abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01);
end

RMSEall = zeros([length(wL) length(wM)]);
markerPlotSpeed = 'sod';
wvMeanAll = zeros([10 3]);
wvPredAll = zeros([10 3]);

for l = 1:length(wL)
    for k = 1:length(wM)
        [~, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCtestWvInFocusMeanZspatFilterPlotHelper(subjNum,defocus875,rgbAll,optDistAll,[wL(l) wM(k) wS]);
        optDistTag = imresize(optDistUnq',size(defocus875mean),'nearest');
        [pFit,RMSE(k)] = ARCfitLagLead(defocus875predTmp(:),defocus875mean(:),optDistTag(:),true,objFunc);
        
        if k==1
            [pFitFlat,RMSEflat(k)] = ARCfitLagLead(optDistTag(:),defocus875mean(:),optDistTag(:),true,objFunc);
            RMSE0 = load([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResults0.mat'],'RMSEall');
            RMSElum = min(RMSE0.RMSEall(1,:));
        end
        
        defocus875pred = [];
        defocus875mean2fit = [];
        figure;
        set(gcf,'Position',[437 314 729 420]);
        subplot(1,2,1);
        hold on;
        for i = 1:size(defocus875predTmp,2)
            hold on;
            dfPred1to5 = -(defocus875predTmp(ind(1:5),i)-optDistUnq(i)*pFit(1)-pFit(2));
            wvPred1to5 = humanWaveDefocusInvertARC(875,-(dfPred1to5+optDistUnq(i)),subjNum);
            dfMean1to5 = -defocus875mean(ind(1:5),i);
            wvMean1to5 = humanWaveDefocusInvertARC(875,-(dfMean1to5+optDistUnq(i)),subjNum);     
            plot(1:5,wvPred1to5,'k-');
            % plot(1:length(ind),defocus875mean(ind,i),'k-');
            for j = 1:5
                plot(j,wvMean1to5(j),['k' markerPlotSpeed(i)],'MarkerFaceColor',conditionsOrderedNorm(j,:), ...
                     'MarkerSize',10);
                defocus875mean2fit(j,i) = defocus875mean(ind(j),i);
            end
            title(['Subject ' num2str(subjNum) ', Distance = ' num2str(optDistUnq(i)) ', Weights = [' num2str(wL(l)) ' ' num2str(wM(k))]);
            title(['RMSE = ' num2str(RMSE,3) ', RMSE_{flat} = ' num2str(RMSEflat,3) ', RMSE_{lum} = ' num2str(RMSElum,3)])
            wvMeanAll(1:5,i) = wvMean1to5;
            wvPredAll(1:5,i) = wvPred1to5;
        end
        set(gca,'FontSize',15);
        set(gca,'XTick',[]);
        xlabel('Condition');
        ylabel('Defocus at 875nm (D)');
        xlim([0 6]);

        subplot(1,2,2);
        for i = 1:size(defocus875predTmp,2)
            hold on;
            dfPred6to10 = -(defocus875predTmp(ind(6:10),i)-optDistUnq(i)*pFit(1)-pFit(2));
            wvPred6to10 = humanWaveDefocusInvertARC(875,-(dfPred6to10+optDistUnq(i)),subjNum);
            dfMean6to10 = -defocus875mean(ind(6:10),i);
            wvMean6to10 = humanWaveDefocusInvertARC(875,-(dfMean6to10+optDistUnq(i)),subjNum);
            dfMean11 = -defocus875mean(ind(11),i);
            wvMean11 = humanWaveDefocusInvertARC(875,-(dfMean11+optDistUnq(i)),subjNum);
            % plot([0 length(ind)],optDistUnq(i)-(optDistUnq(i).*pFitFlat(1)+pFitFlat(2)).*[1 1],'k--','LineWidth',1);
            plot(1:5,wvPred6to10,'k-');
            defocus875pred(:,i) = defocus875predTmp(ind,i)-optDistUnq(i)*pFit(1)-pFit(2);
            % plot(1:length(ind),defocus875mean(ind,i),'k-');
            for j = 6:11
                if j<11
                    plot(j-5,wvMean6to10(j-5),['k' markerPlotSpeed(i)],'MarkerFaceColor',conditionsOrderedNorm(j,:), ...
                         'MarkerSize',10);
                end
                defocus875mean2fit(j,i) = defocus875mean(ind(j),i);
            end
            % plot(6,-(defocus875predTmp(ind(11),i)-optDistUnq(i)*pFit(1)-pFit(2)),'kp','MarkerSize',12);
            wvMeanAll(6:10,i) = wvMean6to10;
            wvPredAll(6:10,i) = wvPred6to10;  
            wvMeanAll(11,i) = wvMean11;
        end
        set(gca,'FontSize',15);
        set(gca,'XTick',[]);
        xlabel('Condition');
        ylabel('Defocus at 875nm (D)');
        title(['RMSE = ' num2str(RMSE,3) ', RMSE_{flat} = ' num2str(RMSEflat,3) ', RMSE_{lum} = ' num2str(RMSElum,3)]);
        xlim([0 6]);
        % saveas(gcf,[coneWeightsFolder 'LplusMminusSpredCont' num2str(subjNum) 'weight' num2str(k)],'png');
        display(['Weights = [' num2str(wL(l)) ' ' num2str(wM(k)) ' ' num2str(wS)]);

    end
end

saveas(gcf,[modelCompFolder 'fitStackSpatFilter' modelType num2str(subjNum)],'png');

errorIndividual = defocus875mean2fit(:)-defocus875pred(:);
for i = 1:200
   [stdTmp(i),LLtmp(i)] = ARCfitStdGauss(errorIndividual);
end
[~,bestInd] = min(LLtmp);
estResidualStd = stdTmp(bestInd);
LL = sum(log(normpdf(defocus875mean2fit(:),defocus875pred(:),estResidualStd)));
aic = 2*nParams-2*LL;

end