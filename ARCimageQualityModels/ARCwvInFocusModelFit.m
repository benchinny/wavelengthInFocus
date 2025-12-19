%% LOADING DATA

function ARCwvInFocusModelFit(subjNum,dataPath)

objFunc = 'RMS'; % OBJECTIVE FUNCTION FOR FITTING
% LM     : LUMINANCE MODEL WITH FITTED WEIGHTS
% LMS    : BLUE-YELLOW MODEL WITH FITTED WEIGHTS
% LminusM: RED-GREEN MODEL WITH FITTED WEIGHTS
LMorLMSorLminusM = 'LMS'; 

% LIST OF ALL SUBJECTS
subjNumListAll = [1 3 5 10 16 17 18 20];
% FIND subjNum POSITION IN ARRAY
indLCA = find(subjNumListAll==subjNum);
% LOAD PRE-SAVED LCA PARAMETERS
load(fullfile(dataPath,'data','PresavedFigureData','LCAparams.mat'),'q1bestAll','q2bestAll','q3bestAll');
q1 = q1bestAll(indLCA);
q2 = q2bestAll(indLCA);
q3 = q3bestAll(indLCA);

% TAG SUBJECT FILE / BLOCK NUMBERS
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
        [ZernikeTable, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumAll(k),trialNumAll(l),dataPath);

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

% EXCLUDE DATA FOR WHICH PARTICIPANT WAS ACCOMMODATING OUTSIDE OF
% VISIBLE RANGE
% DEFOCUS AT 550NM
defocus550 = defocus875+humanWaveDefocusParameterizedARC(550,875,q1,q2,q3);
% DEVIATION FROM STIMULUS DISTANCE AT BOTH 550NM AND 875NM
diffFromOptDist875 = defocus875-optDistAll;
diffFromOptDist550 = defocus550-optDistAll;
% 'GOOD INDICES' AT WHICH SUBJECT IS ACCOMMODATING WITHIN VISIBLE RANGE
indGood = abs(diffFromOptDist550)<2 & ...
          humanWaveDefocusInvertParameterizedARC(875,diffFromOptDist875,q1,q2,q3)>380 & ...
          humanWaveDefocusInvertParameterizedARC(875,diffFromOptDist875,q1,q2,q3)<780;
defocus875 = defocus875(indGood);
rgbAll = rgbAll(indGood,:);
optDistAll = optDistAll(indGood);

%% SEARCH INDIVIDUAL CONE WEIGHTS

if strcmp(LMorLMSorLminusM,'LMS')
    % WEIGHT VALUES FOR GRID SEARCH
    wLM = 0.4:0.05:1.4; % ratio of (L+M) to S
    wLprop = 0.25:(0.1/3):0.85; % ratio of L to M
    if subjNum==5
        wS = 0.5;
    elseif subjNum==20
        wS = 0.25;
    else
        wS = 1;
    end
    modelResultsFilename = 'wvInFocusModelResultsDonutx2';
elseif strcmp(LMorLMSorLminusM,'LM')
    % WEIGHT VALUES FOR GRID SEARCH
    wLM = [0.5 1]; % ratio of (L+M) to S
    wLprop = 0.25:(0.1/3):0.85; % ratio of L to M    
    wS = 0;
    modelResultsFilename = 'wvInFocusModelResults';
elseif strcmp(LMorLMSorLminusM,'LminusM')
    % WEIGHT VALUES FOR GRID SEARCH
    wLM = [0.5 1]; % ratio of (L+M) to S
    wLprop = 0.25:(0.1/3):0.85; % ratio of L to M    
    wS = 0;
    modelResultsFilename = 'wvInFocusModelResultsLminusM';
else
    error('Specify valid model type!');
end

coneWeightsFolder = fullfile(dataPath,'data','coneWeightsErrorSpatFilter','colorMechPredictions');

RMSEall = zeros([length(wLM) length(wLprop)]); % INITIALIZE ERROR SURFACE

for l = 1:length(wLM) % LOOP OVER RATIO OF L+M TO S
    parfor k = 1:length(wLprop) % LOOP OVER L TO M RATIO
        % CONVERTING TO WEIGHTS ON L AND M
        if strcmp(LMorLMSorLminusM,'LMS') || strcmp(LMorLMSorLminusM,'LM')
            wL = wLM(l)*wLprop(k);
            wM = wLM(l)-wL;
        elseif strcmp(LMorLMSorLminusM,'LminusM')
            wL = wLM(l)*wLprop(k);
            wM = -(wLM(l)-wL); 
        end
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

save([coneWeightsFolder 'S' num2str(subjNum) modelResultsFilename num2str(round(-wS*10)) '.mat'],'RMSEall','wS','wLM','wLprop');

end
