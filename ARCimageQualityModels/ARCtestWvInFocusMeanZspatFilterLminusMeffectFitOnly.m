%% LOADING DATA

function ARCtestWvInFocusMeanZspatFilterLminusMeffectFitOnly(subjNum,dataPath)

objFunc = 'RMS'; % ERROR FUNCTION FOR FITTING

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
optDistAll = []; % OPTICAL DISTANCES OF STIMULI
rgbAll = []; % STIMULUS COLOR CONDITIONS

% LOAD PARTICIPANTS' DATA TO FIT
for k = 1:length(blockNumAll) % LOOP OVER BLOCK NUMBER
    AFCp = ARCloadFileBVAMS(subjNum+10,blockNumAll(k),dataPath); % LOAD EXPERIMENT FILE
    optDistAll = [optDistAll; AFCp.meanv00./1.2255]; % STACK OPTICAL DISTANCES (1.2255 IS SCALE FACTOR FOR BVAMS)
    rgbAll = [rgbAll; AFCp.rgb100]; % STACK COLOR CONDITIONS
    for l = 1:length(trialNumAll) % LOOP OVER TRIALS
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumAll(k),trialNumAll(l),0,dataPath);

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS = struct;
        indBadPupil = table2array(ZernikeTable(:,5))==0; % FIND BLINKS
        PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0; % FIND BLINKS IN DEFOCUS VECTOR
        meanC = mean(c(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS
        % STANDARD SCALAR FOR CONVERTING ZERNIKE COEFFICIENT TO EQUIVALENT
        % DEFOCUS
        defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
        defocus875(end+1,:) = meanC(4)./defocusCorrectionFactor;
    end
end

%% SEARCH INDIVIDUAL CONE WEIGHTS

% WEIGHTS FOR GRID SEARCH
wLM = [0.5 1]; % TOTAL WEIGHT ON L AND M (SHOULDN'T AFFECT FIT)
wLprop = 0:0.025:1; % L TO M RATIO

if ispc
    slash = '\';
else
    slash = '/';
end
coneWeightsFolder = [dataPath 'data' slash 'coneWeightsErrorSpatFilter' slash 'colorMechPredictions' slash];

RMSEall = zeros([length(wLM) length(wLprop)]); % INITIALIZE ERROR SURFACE MATRIX

for l = 1:length(wLM) % LOOP OVER TOTAL L AND M WEIGHT
    parfor k = 1:length(wLprop) % LOOP OVER L PROPORTION
        % CONVERT TO INDIVIDUAL WEIGHTS
        wL = wLM(l)*wLprop(k);
        wM = -(wLM(l)-wL);
        wS = 0;
        % USE HELPER FUNCTION TO GENERATE PREDICTIONS OF DEFOCUS
        [~, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCtestWvInFocusMeanZspatFilterPlotHelper(subjNum,defocus875,rgbAll,optDistAll,[wL wM wS],dataPath);
        % TAG OPTICAL DISTANCES FOR FITTING LAG AND LEAD PARAMETERS
        optDistTag = imresize(optDistUnq',size(defocus875mean),'nearest');
        % FIT LAGS AND LEADS FREE PARAMETERS
        [pFit,RMSE(k)] = ARCfitLagLead(defocus875predTmp(:),defocus875mean(:),optDistTag(:),true,objFunc);
        
        display(['Weights = [' num2str(wL) ' ' num2str(wM) ' ' num2str(wS)]);
    end
    RMSEall(l,:) = RMSE;
end

save([coneWeightsFolder 'S' num2str(subjNum) 'wvInFocusModelResultsLminusM.mat'],'RMSEall','wLM','wLprop');

end
