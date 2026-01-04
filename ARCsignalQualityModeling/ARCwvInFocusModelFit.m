%% LOADING DATA

function ARCwvInFocusModelFit(subjNum,modelType,sigQualType,dataPath)

% subjNum: subject number. Valid numbers: 1, 3, 5, 10, 16, 17, 18, 20
% modelType options
%                   LM     : LUMINANCE MODEL WITH FITTED WEIGHTS
%                   LMS    : BLUE-YELLOW MODEL WITH FITTED WEIGHTS
%                   LminusM: RED-GREEN MODEL WITH FITTED WEIGHTS
% sigQualType options
%                   'xcorr' : cross-correlation metric (main paper)
%                   'strehl': normalized Strehl (supplement)
%                   'deltapass': Finch et al. model, which is bandpass at a
%                                single frequency (hence 'delta pass')
% dataPath: local directory for data

rng(1); % INITIALIZE SAME RANDOM SEED

objFunc = 'RMS'; % OBJECTIVE FUNCTION FOR FITTING

% LOAD DEFOCUS VALUES, COLOR CONDITIONS, AND OPTICAL DISTANCES
[defocus875,rgbAll,optDistAll,~,~] = ARCnlzLoadDefocusAbb(subjNum,dataPath);

%% SEARCH INDIVIDUAL CONE WEIGHTS

if strcmp(modelType,'LMS')
    % WEIGHT VALUES FOR GRID SEARCH
    wLM = 0.4:0.05:1.4; % ratio of (L+M) to S
    wLprop = 0.25:(0.1/3):0.85; % ratio of L to M
    % FOR THESE TWO SUBJECTS, IF wS = 1, THE SEARCH WILL BUMP UP AGAINST
    % THE EDGE OF THE GRID
    if subjNum==5
        wS = -0.5;
    elseif subjNum==20
        wS = -0.25;
    else
        wS = -1;
    end
    modelResultsFilename = 'wvInFocusModelResultsLMS';
elseif strcmp(modelType,'LM')
    % WEIGHT VALUES FOR GRID SEARCH
    wLM = [0.5 1]; % ratio of (L+M) to S
    wLprop = 0.25:(0.1/3):0.85; % ratio of L to M    
    wS = 0;
    modelResultsFilename = 'wvInFocusModelResults';
elseif strcmp(modelType,'LminusM')
    % WEIGHT VALUES FOR GRID SEARCH
    wLM = [0.5 1]; % ratio of (L+M) to S
    wLprop = 0.25:(0.1/3):0.85; % ratio of L to M    
    wS = 0;
    modelResultsFilename = 'wvInFocusModelResultsLminusM';
elseif strcmp(modelType,'Lum')
    wLM = 1;
    wLprop = 0.72;
    wS = 0;
    modelResultsFilename = 'wvInFocusModelResultsLum';
else
    error('Specify valid model type!');
end

% SPECIFY STRING TO APPEND TO FILE NAME
if strcmp(sigQualType,'xcorr')
    metricName = '';
elseif strcmp(sigQualType,'strehl')
    metricName = 'strehl';
elseif strcmp(sigQualType,'deltapass')
    metricName = 'deltapass';    
else
    error('Specify valid string for sigQualType: either xcorr, strehl, or deltapass');
end

% PATH TO CONE WEIGHTS
coneWeightsFolder = fullfile(dataPath,'data','coneWeights');

RMSEall = zeros([length(wLM) length(wLprop)]); % INITIALIZE ERROR SURFACE
pFitAll = zeros([length(wLM) length(wLprop) 2]);

for l = 1:length(wLM) % LOOP OVER RATIO OF L+M TO S
    parfor k = 1:length(wLprop) % LOOP OVER L TO M RATIO
        % CONVERTING TO WEIGHTS ON L AND M
        if strcmp(modelType,'LMS') || strcmp(modelType,'LM') || strcmp(modelType,'Lum')
            wL = wLM(l)*wLprop(k);
            wM = wLM(l)-wL;
        elseif strcmp(modelType,'LminusM')
            wL = wLM(l)*wLprop(k);
            wM = -(wLM(l)-wL);     
        end
        % GENERATE PREDICTIONS OF DEFOCUS USING HELPER FUNCTION
        [~, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCwvInFocusModelHelper(subjNum,defocus875,rgbAll,optDistAll,[wL wM wS],sigQualType,dataPath);
        % TAG EVERY TRIAL BY OPTICAL DISTANCE FOR FITTING LAGS AND LEADS
        optDistTag = imresize(optDistUnq',size(defocus875mean),'nearest');
        % FIT LAGS AND LEADS
        [pFit(k,:),RMSE(k)] = ARCfitLagLead(defocus875predTmp(:),defocus875mean(:),optDistTag(:),true,objFunc);

        display(['Weights = [' num2str(wL) ' ' num2str(wM) ' ' num2str(wS)]);
    end
    RMSEall(l,:) = RMSE;
    pFitAll(l,:,:) = pFit;
end

save(fullfile(coneWeightsFolder,['S' num2str(subjNum) modelResultsFilename metricName num2str(round(-wS*10)) '.mat']),'RMSEall','wS','wLM','wLprop','pFitAll');

end
