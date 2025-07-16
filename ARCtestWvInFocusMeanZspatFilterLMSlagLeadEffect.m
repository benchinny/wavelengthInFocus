function [waveInFocus3raw, waveInFocus3, wLMS, lagLead] = ARCtestWvInFocusMeanZspatFilterLMSlagLeadEffect(subjNum,modelType,stimNum)

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

for i = 1:size(conditionsOrderedNorm,1)
    ind(i) = find(abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
                  abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
                  abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01);
end

RMSEall = zeros([length(wL) length(wM)]);
markerPlotSpeed = 'sod';

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
        
        lagLead = [];
        waveInFocus3 = [];

        for i = 1:size(defocus875predTmp,2)
            waveInFocus3(i) = humanWaveDefocusInvertARC(875,-(-(defocus875predTmp(ind(stimNum),i)-optDistUnq(i)*pFit(1)-pFit(2))+optDistUnq(i)),subjNum);
            lagLead(i) = optDistUnq(i)*pFit(1)-pFit(2);
        end
        waveInFocus3raw = humanWaveDefocusInvertARC(875,-(-(defocus875predTmp(ind(stimNum),i))+optDistUnq(i)),subjNum);
        wLMS = [wL(l) wM(k) wS];
        display(['Weights = [' num2str(wL(l)) ' ' num2str(wM(k)) ' ' num2str(wS)]);
    end
end

end