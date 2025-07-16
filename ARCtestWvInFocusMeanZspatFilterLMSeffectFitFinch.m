%% LOADING DATA

function ARCtestWvInFocusMeanZspatFilterLMSeffectFitFinch(wS)

subjNum = 21;
objFunc = 'RMS';

%% SEARCH INDIVIDUAL CONE WEIGHTS

% wLM = 0.4:0.05:1.4;
% wLprop = 0.25:(0.1/3):0.85;
% wLM = 1.4;
% wLprop = 0.25:(0.1/3):0.85;

if subjNum==1
    wLM = 0.8;
    wLprop = 0.6833;
elseif subjNum==10
    wLM = 1.35;
    wLprop = 0.85;
elseif subjNum==16
    wLM = 0.8;
    wLprop = 0.2833;
elseif subjNum==21
    % wLM = 1./[1.00 0.90 0.80 0.70 0.60 0.50 0.40 0.30 0.20 0.10];
    % wLprop = 0.25:0.05:0.85; 
    wLM = 1;
    wLprop = 0.65;
end

coneWeightsFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsErrorSpatFilter/colorMechPredictions/';
finchDataFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/';
colorCell = {'RedVsGreen' 'RedVsBlue' 'RedVsViolet' 'OrangeVsBlue' 'OrangeVsViolet' 'GreenVsViolet'};
load([finchDataFolder 'Finch_data_all.mat'],'FinchData');

defocus550 = FinchData.Accommodation;
optDistUnq = 0;

RMSEall = zeros([length(wLM) length(wLprop)]);

for l = 1:length(wLM)
    for k = 1:length(wLprop)
        wL = wLM(l)*wLprop(k);
        wM = wLM(l)-wL;
        [~, defocus550mean, defocus550predTmp, ~, ~] = ARCtestWvInFocusMeanZspatFilterFinchPlotHelper(subjNum,flipud(defocus550),flipud(FinchData.LumRatio),optDistUnq,[wL wM wS]);
        optDistTag = imresize(optDistUnq',size(defocus550mean),'nearest');
        [pFit,RMSE(k)] = ARCfitLagLead(defocus550predTmp(:),defocus550mean(:),optDistTag(:),true,objFunc);
        if k==1
            [pFitFlat,RMSEflat(k)] = ARCfitLagLead(optDistTag(:),defocus550mean(:),optDistTag(:),true,objFunc);
        end
        
        display(['Weights = [' num2str(wL) ' ' num2str(wM) ' ' num2str(wS)]);
    end
    RMSEall(l,:) = RMSE;
end

% save([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResultsFinch' num2str(round(-wS*10)) '.mat'],'RMSEall','wS','wLM','wLprop');

end
