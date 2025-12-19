function [aic, pFit, wvMeanAll, wvPredAll, dfPredPurple, wLMmin, wLpropMin] = ARCwvInFocusModelSort(subjNum,modelType,dataPath)

% WHERE THE PRE-FIT CONE WEIGHTS ARE STORED
coneWeightsFolder = fullfile(dataPath,'data','coneWeightsErrorSpatFilter','colorMechPredictions');

objFunc = 'RMS'; % OBJECTIVE FUNCTION FOR EVALUATING FIT

% LIST OF ALL SUBJECTS
subjNumListAll = [1 3 5 10 16 17 18 20];
% FIND subjNum POSITION IN ARRAY
indLCA = find(subjNumListAll==subjNum);
% LOAD PRE-SAVED LCA PARAMETERS
load(fullfile(dataPath,'data','PresavedFigureData','LCAparams.mat'),'q1bestAll','q2bestAll','q3bestAll');
q1 = q1bestAll(indLCA);
q2 = q2bestAll(indLCA);
q3 = q3bestAll(indLCA);

% LOAD CONE WEIGHTS
if strcmp(modelType,'LMS') % IF BLUE-YELLOW OPPONENT MODEL
    % THIS VALUE WAS GOOD TO SEARCH OVER A WIDE RANGE OF (W_L+W_M)/W_S
    % VALUES
    wS = -1; 
    % FOR THE FOLLOWING TWO SUBJECTS, A MUCH HIGHER (W_L+W_M)/W_S RATIO WAS
    % REQUIRED, SO wS WAS FIXED TO A LOWER VALUE SO THAT THE GRID SEARCH
    % COULD SEARCH OVER HIGHER RATIOS
    if subjNum==20
        wS = -0.25;
    end
    if subjNum==5
        wS = -0.5;
    end
    % LOAD PRE-FIT CONE WEIGHTS
    load(fullfile(coneWeightsFolder,['S' num2str(subjNum) 'wvInFocusModelResultsDonutx2' num2str(round(-wS*10)) '.mat']),'RMSEall','wS','wLM','wLprop');
    
    % MAKE A MESHGRID FOR EASILY FINDING THE BEST FIT PRE-GENERATED PARAMETERS
    [wLpropGrid,wLMgrid] = meshgrid(wLprop,wLM);
    % INDEX IN GRID FOR BEST PARAMETERS
    indMin = RMSEall==min(RMSEall(:));
    % BEST FIT PARAMETERS
    wLMmin = wLMgrid(indMin);
    wLpropMin = wLpropGrid(indMin);
    % CONVERT BEST FIT PARAMETERS TO L, M, AND S WEIGHTS
    wL = wLMmin*wLpropMin;
    wM = wLMmin-wL;
    nParams = 4; % FOR CALCULATING AIC LATER
end

if strcmp(modelType,'LM') % IF LUMINANCE MODEL WITH FREE WEIGHTS
    wS = 0; % ALWAYS 0 BY DEFINITION
    % LOAD PRE-FIT CONE WEIGHTS
    load(fullfile(coneWeightsFolder,['S' num2str(subjNum) 'wvInFocusModelResults' num2str(round(-wS*10)) '.mat']),'RMSEall','wS','wLM','wLprop');
    
    % MAKE A MESHGRID FOR EASILY FINDING THE BEST FIT PRE-GENERATED PARAMETERS
    [wLpropGrid,wLMgrid] = meshgrid(wLprop,wLM);
    % INDEX IN GRID FOR BEST PARAMETERS
    indMin = RMSEall==min(RMSEall(:));
    % BEST FIT PARAMETERS
    wLMmin = wLMgrid(indMin);
    wLpropMin = wLpropGrid(indMin);
    % CONVERT BEST FIT PARAMETERS TO L AND M WEIGHTS
    wL = wLMmin*wLpropMin;
    wM = wLMmin-wL;
    nParams = 3; % FOR CALCULATING AIC LATER
end

if strcmp(modelType,'LminusM') % IF RED-GREEN OPPONENT MODEL
    wS = 0; % ALWAYS 0 BY DEFINITION
    % LOAD PRE-FIT WEIGHTS
    load(fullfile(coneWeightsFolder,['S' num2str(subjNum) 'wvInFocusModelResultsLminusM.mat']),'RMSEall','wLM','wLprop');
    
    % MAKE A MESHGRID FOR EASILY FINDING THE BEST FIT PRE-GENERATED PARAMETERS
    [wLpropGrid,wLMgrid] = meshgrid(wLprop,wLM);
    % INDEX IN GRID FOR BEST PARAMETERS
    indMin = RMSEall==min(RMSEall(:));
    % BEST FIT PARAMETERS
    wLMmin = wLMgrid(indMin);
    wLpropMin = wLpropGrid(indMin);
    % CONVERT TO L AND M WEIGHTS
    wL = wLMmin*wLpropMin;
    wM = -(wLMmin-wL);    
    nParams = 3; % FOR CALCULATING AIC LATER
end

% WE DON'T REALLY USE THIS IN THE PAPER
if strcmp(modelType,'Lum') % IF LUMINANCE MODEL WITH V-LAMBDA
    wL = 0.72;
    wM = 0.28;
    wS = 0;
    nParams = 2;
end

% SPECIFY BLOCKS FOR DIFFERENT SUBJECTS
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

trialNumAll = 1:36; % SAME NUMBER OF TRIALS FOR ALL SUBJECTS

defocus875 = []; % DEFOCUS ABERRATION AT 875NM
optDistAll = []; % OPTICAL DISTANCES
rgbAll = []; % COLOR CONDITIONS

% LOADING EMPIRCAL DATA TO PLOT WITH PREDICTIONS
for k = 1:length(blockNumAll) % LOOP OVER BLOCKS
    % LOAD EXPERIMENT DATA
    AFCp = ARCloadFileBVAMS(subjNum+10,blockNumAll(k),dataPath);
    % 1.2255 IS THE SCALE FACTOR FOR CONVERTING BVAMS POWER CHANGE TO
    % DEFOCUS AT EYE
    optDistAll = [optDistAll; AFCp.meanv00./1.2255];
    rgbAll = [rgbAll; AFCp.rgb100]; % CONCATENATE COLOR CONDITIONS
    for l = 1:length(trialNumAll) % LOOP OVER TRIALS
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumAll(k),trialNumAll(l),dataPath);

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS = struct;
        indBadPupil = table2array(ZernikeTable(:,5))==0; % FIND BLINKS
        PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0; % FIND BLINKS IN DEFOCUS VECTOR
        c(indBad,4) = mean(c(~indBad,4)); % REPLACE BLINK TIME POINTS
        meanC = mean(c(1:end,:),1); % TAKE MEAN OF COEFFICIENTS
        % STANDARD CORRECTION FACTOR FOR GOING FROM COEFFICIENT TO
        % EQUIVALENT DEFOCUS
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

rgbUnq = unique(rgbAll,'rows'); % UNIQUE RGB VALUES

% NORMALIZE RGB VALUES SO MAX LUMINANCE IS 1
rgbLumNorm = [];
rgbLumNorm(:,1) = (rgbUnq(:,1).^2.5)./0.2442;
rgbLumNorm(:,2) = (rgbUnq(:,2).^2.7)./0.1037;
rgbLumNorm(:,3) = (rgbUnq(:,3).^2.3)./1;
rgbLumNorm(rgbLumNorm>1) = 1;

% ORDER OF CONDITIONS FOR PLOTTING: MORE BLUE TO MORE RED WITHOUT GREEN,
% THEN MORE BLUE TO MORE RED WITH GREEN
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

% INDICES FOR SORTING CONDITIONS ACCORDING TO PLOT
for i = 1:size(conditionsOrderedNorm,1)
    ind(i) = find(abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
                  abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
                  abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01);
end

wvMeanAll = zeros([10 3]);
wvPredAll = zeros([10 3]);

for l = 1:length(wL)
    for k = 1:length(wM)
        % IMPORTANT THINGS HAPPENING IN HELPER FUNCTION: GENERATE
        % PREDICTIONS OF DEFOCUS FOR EACH CONDITION
        [~, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCtestWvInFocusMeanZspatFilterPlotHelper(subjNum,defocus875,rgbAll,optDistAll,[wL(l) wM(k) wS],dataPath);
        % VECTOR OF OPTICAL DISTANCES FOR TAGGING FOR FITTING LAGS AND
        % LEADS
        optDistTag = imresize(optDistUnq',size(defocus875mean),'nearest');
        % FIT LAGS AND LEADS
        [pFit,RMSE(k)] = ARCfitLagLead(defocus875predTmp(:),defocus875mean(:),optDistTag(:),true,objFunc);
        
        % SORTING DATA FOR PLOTTING
        defocus875pred = [];
        defocus875mean2fit = [];
        % 'NO GREEN' CONDITIONS ARE INDICES 1 TO 5
        for i = 1:size(defocus875predTmp,2)
            % PREDICTIONS (NEED TO CONVERT FROM DEFOCUS TO WAVELENGTH)
            dfPred1to5 = -(defocus875predTmp(ind(1:5),i)-optDistUnq(i)*pFit(1)-pFit(2));
            wvPred1to5 = humanWaveDefocusInvertParameterizedARC(875,-(dfPred1to5+optDistUnq(i)),q1,q2,q3);
            % ACTUAL (NEED TO CONVERT FROM DEFOCUS TO WAVELENGTH)
            dfMean1to5 = -defocus875mean(ind(1:5),i);
            wvMean1to5 = humanWaveDefocusInvertParameterizedARC(875,-(dfMean1to5+optDistUnq(i)),q1,q2,q3);
            if i==2 % STORE AND RETURN PREDICTION FOR USE IN ACUITY MODELING
               dfPredPurple = dfPred1to5(3);
            end
            for j = 1:5
                % USE THIS VALUE TO CALCULATE AIC LATER
                defocus875mean2fit(j,i) = defocus875mean(ind(j),i);
            end
            % STORE MEAN DATA AND PREDICTIONS
            wvMeanAll(1:5,i) = wvMean1to5;
            wvPredAll(1:5,i) = wvPred1to5;
        end
        % 'GREEN CONDITIONS ARE INDICES 6 TO 10
        for i = 1:size(defocus875predTmp,2)
            % PREDICTIONS (NEED TO CONVERT FROM DEFOCUS TO WAVELENGTH)
            dfPred6to10 = -(defocus875predTmp(ind(6:10),i)-optDistUnq(i)*pFit(1)-pFit(2));
            wvPred6to10 = humanWaveDefocusInvertParameterizedARC(875,-(dfPred6to10+optDistUnq(i)),q1,q2,q3);
            % ACTUAL (NEED TO CONVERT FROM DEFOCUS TO WAVELENGTH)
            dfMean6to10 = -defocus875mean(ind(6:10),i);
            wvMean6to10 = humanWaveDefocusInvertParameterizedARC(875,-(dfMean6to10+optDistUnq(i)),q1,q2,q3);
            % BREAK OUT CONTROL CONDITION TO RETURN SEPARATELY
            dfMean11 = -defocus875mean(ind(11),i);
            wvMean11 = humanWaveDefocusInvertParameterizedARC(875,-(dfMean11+optDistUnq(i)),q1,q2,q3);
            % STORE ALL PREDICTIONS (ONLY NEED TO DO THIS IN SECOND LOOP)
            defocus875pred(:,i) = defocus875predTmp(ind,i)-optDistUnq(i)*pFit(1)-pFit(2);
            % USE THIS VALUE TO CALCULATE AIC LATER
            for j = 6:11
                defocus875mean2fit(j,i) = defocus875mean(ind(j),i);
            end
            % STORE MEAN DATA AND PREDICTIONS
            wvMeanAll(6:10,i) = wvMean6to10;
            wvPredAll(6:10,i) = wvPred6to10;  
            wvMeanAll(11,i) = wvMean11;
        end
        display(['Weights = [' num2str(wL(l)) ' ' num2str(wM(k)) ' ' num2str(wS)]);
    end
end

% CALCULATING AIC--FIRST CALCULATE DIFFERENCE BETWEEN MEAN DEFOCUS AND
% PREDICTIONS
errorIndividual = defocus875mean2fit(:)-defocus875pred(:);
% ESTIMATE STD DEV OF RESIDUALS
for i = 1:200
   [stdTmp(i),LLtmp(i)] = ARCfitStdGauss(errorIndividual);
end
[~,bestInd] = min(LLtmp);
estResidualStd = stdTmp(bestInd);
% CALCULATE LOG-LIKELIHOOD, THEN AIC
LL = sum(log(normpdf(defocus875mean2fit(:),defocus875pred(:),estResidualStd)));
aic = 2*nParams-2*LL;

end