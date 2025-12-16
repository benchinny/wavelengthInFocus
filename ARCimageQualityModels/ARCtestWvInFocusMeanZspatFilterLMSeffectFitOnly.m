%% LOADING DATA

function ARCtestWvInFocusMeanZspatFilterLMSeffectFitOnly(subjNum,wS,dataPath)

objFunc = 'RMS'; % OBJECTIVE FUNCTION FOR FITTING

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

trialNumAll = 1:36; % ALL SUBJECTS HAVE SAME NUMBER OF TRIALS

defocus875 = []; % DEFOCUS AT 875NM
optDistAll = []; % STIMULUS DISTANCES
rgbAll = []; % COLOR CONDITIONS

% LOAD DATA TO FIT
for k = 1:length(blockNumAll) % LOOP OVER BLOCKS
    AFCp = ARCloadFileBVAMS(subjNum+10,blockNumAll(k),dataPath); % LOAD CONDITION FILE
    optDistAll = [optDistAll; AFCp.meanv00./1.2255]; % CONCATENATE OPTICAL DISTANCES (1.2255 SCALE FACTOR)
    rgbAll = [rgbAll; AFCp.rgb100]; % CONCATENATE COLOR CONDITIONS
    for l = 1:length(trialNumAll) % FOR EACH TRIAL
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumAll(k),trialNumAll(l),0,dataPath);

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS = struct;
        indBadPupil = table2array(ZernikeTable(:,5))==0; % IDENTIFY BLINKS
        PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        % ALL ZERNIKE COEFFICIENTS
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0; % FIND BLINKS IN DEFOCUS VECTOR
        meanC = mean(c(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS WITHOUT BLINKS
        % STANDARD CORRECTION TO CONVERT TO EQUIVALENT DEFOCUS
        defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
        defocus875(end+1,:) = meanC(4)./defocusCorrectionFactor;
    end
end

%% SEARCH INDIVIDUAL CONE WEIGHTS

% WEIGHT VALUES FOR GRID SEARCH
wLM = 0.4:0.05:1.4; % ratio of (L+M) to S
wLprop = 0.25:(0.1/3):0.85; % ratio of L to M

if ispc
    slash = '\';
else
    slash = '/';
end
coneWeightsFolder = [dataPath 'data' slash 'coneWeightsErrorSpatFilter' slash 'colorMechPredictions' slash];

RMSEall = zeros([length(wLM) length(wLprop)]); % INITIALIZE ERROR SURFACE

for l = 1:length(wLM) % LOOP OVER RATIO OF L+M TO S
    parfor k = 1:length(wLprop) % LOOP OVER L TO M RATIO
        wL = wLM(l)*wLprop(k);
        wM = wLM(l)-wL;
        % GENERATE PREDICTIONS OF DEFOCUS USING HELPER FUNCTION
        [~, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCtestWvInFocusMeanZspatFilterPlotHelper(subjNum,defocus875,rgbAll,optDistAll,[wL wM wS],dataPath);
        % TAG EVERY TRIAL BY OPTICAL DISTANCE FOR FITTING LAGS AND LEADS
        optDistTag = imresize(optDistUnq',size(defocus875mean),'nearest');
        % FIT LAGS AND LEADS
        [pFit,RMSE(k)] = ARCfitLagLead(defocus875predTmp(:),defocus875mean(:),optDistTag(:),true,objFunc);
        
        display(['Weights = [' num2str(wL) ' ' num2str(wM) ' ' num2str(wS)]);
    end
    RMSEall(l,:) = RMSE;
end

save([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResultsDonutx2' num2str(round(-wS*10)) '.mat'],'RMSEall','wS','wLM','wLprop');

end
