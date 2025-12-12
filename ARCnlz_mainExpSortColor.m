%% LOAD MAIN EXPERIMENT FILES

function [wvInFocusCell, optDistCnd, rgbLumNormCnd] = ARCnlz_mainExpSortColor(subjNum,dataPath)

% subjNum values for participants who passed screening: S1, S3, S5, S10, S16, S17, S18, S20. 

% This function grabs all raw wavefront data from Experiment 1 for a
% particular subject and returns mean defocus values (both at 550nm and 875nm),
% wavelength-in-focus values, and the optical distances and color
% conditions corresponding to those values. The 'c' variable contains the
% raw traces of every coefficient on each of the 65 Zernike terms. 

% LIST OF ALL POSSIBLE SUBJECTS
subjNumAll = [1 3 5 10 16 17 18 20];
% COMPUTE LCA PARAMETERS FOR ALL SUBJECTS
[q1bestAll, q2bestAll, q3bestAll] = ARCacuAnalysisLCAall(dataPath,0);

if subjNum==1
   blockNums = 11:16;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};   
   subjName = ['S' num2str(subjNum+10) '-OD'];
elseif subjNum==3
   blockNums = 12:17;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};   
   subjName = ['S' num2str(subjNum+10) '-OD'];   
elseif subjNum==5
   blockNums = 3:8;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum+10) '-OD'];   
elseif subjNum==10
   blockNums = 3:8;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum+10) '-OD'];         
elseif subjNum==16
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum+10) '-OD'];      
elseif subjNum==17
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum+10) '-OD'];      
elseif subjNum==18
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum+10) '-OD'];   
elseif subjNum==20
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum+10) '-OD'];      
end

meanC = []; % MEAN OF Z COEFFICIENTS
rgb1all = []; % RGB VALUES
meanv00all = []; % OPTOTUNE VALUE (TO GET STIMULUS DISTANCE, DIVIDE BY 1.2255)
nIndBadTracker = []; % TRACKING BLINKS
c4all = {}; % FOR STORING ALL DEFOCUS TERM VALUES

% GRABBING ZERNIKE VALUES FROM DATA REPOSITORY
for i = 1:length(blockNums)
    blockNumTmp = blockNums(i); % BLOCK REFERS TO BLOCK OF TRIALS
    trialNumsTmp = trialNums{i}; % TRIAL NUMBERS
    AFCp = ARCloadFileBVAMS(subjNum+10,blockNumTmp,dataPath); % LOAD EXPERIMENT FILE
    for j = 1:length(trialNumsTmp) % LOOP OVER TRIALS
        % THIS LOADS THE WAVEFRONT DATA
        [ZernikeTable, ~, ~, ~] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumsTmp(j),0,dataPath);
        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65.
        indBadPupil = table2array(ZernikeTable(:,5))==0; % GET RID OF BLINKS IN PUPIL SIZE VECTOR!
        PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0; % VALUE OF EXACTLY 0 MEANS BLINK
        nIndBadTracker(end+1) = sum(indBad); % TRACK BLINKS
        c(indBad,4) = mean(c(~indBad,4)); % REPLACE BLINKS WITH MEAN
        meanC(end+1,:) = mean(c(1:end,:),1); % TAKE MEAN OF COEFFICIENTS    
        c4all{end+1} = c(:,4); % COEFFICIENTS ON DEFOCUS TERM
    end
    rgb1all = [rgb1all; AFCp.rgb100(trialNumsTmp,:)]; % STORE COLOR CONDITIONS
    meanv00all = [meanv00all; AFCp.meanv00(trialNumsTmp)./1.2255]; % STORE STIM DISTANCES
end

% STANDARD CONVERSION FROM 4TH ZERNIKE COEFFICIENT TO EQUIVALENT DEFOCUS 
defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
% GETTING DEFOCUS AT 550NM and 875NM (875NM IS WHAT WAS MEASURED)
defocusAt550 = humanWaveDefocusParameterizedARC(550,875,q1bestAll(subjNumAll==subjNum),q2bestAll(subjNumAll==subjNum),q3bestAll(subjNumAll==subjNum))+meanC(:,4)./defocusCorrectionFactor;

% SORTING CONDITIONS BY COLOR
lumScaleRGB = [4.0888 9.6669 1]; % SCALE FACTOR FOR NORMALIZING RGB TO 1

gammaRGB = [2.5 2.7 2.3]; % GAMMA EXPONENTS FOR DISPLAY

% NORMALIZE RGB VALUES SO MAXIMUM OF ANY INDIVIDUAL DISPLAY PRIMARY IS 1
rgbLumNorm = [lumScaleRGB(1).*rgb1all(:,1).^gammaRGB(1) lumScaleRGB(2).*rgb1all(:,2).^gammaRGB(2) lumScaleRGB(3).*rgb1all(:,3).^gammaRGB(3)];

% EACH TRIPLET INDICATES THE LUMINANCE OF THE RED, GREEN, AND BLUE PIXELS
% RESPECTIVELY AS A PROPORTION OF THE MAXIMUM VALUE ALLOWED IN THE
% EXPERIMENT
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

optDistToCheckAll = [1.5 2.5 3.5]; % STIMULUS DISTANCES
optDistCnd = []; % FOR SORTING OPTICAL DISTANCES
rgbLumNormCnd = []; % FOR SORTING RGB CONDITION VALUES
wvInFocusCell = {}; % FOR STORING WAVELENGTH IN FOCUS VALUES

for j = 1:length(optDistToCheckAll) % LOOP OVER STIM DISTANCE
    optDistToCheck = optDistToCheckAll(j);
    for i = 1:size(conditionsOrderedNorm,1) % LOOP OVER RGB CONDITIONS
        % GET ALL TRIALS FOR EACH COMBINATION OF COLOR AND DISTANCE
        ind = abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
              abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
              abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01 & ...
              abs(meanv00all-optDistToCheck)<0.01;
        % FOR OUTLIER EXCLUSION
        defocusAt550tmp = defocusAt550(ind);
        % DEFOCUS RELATIVE TO OPTICAL DISTANCE
        diffFromOptDist = defocusAt550tmp-meanv00all(ind);
        % EXCLUDE DATA FOR WHICH PARTICIPANT WAS ACCOMMODATING OUTSIDE OF
        % VISIBLE RANGE
        indGood = abs(diffFromOptDist)<2 & ...
                  humanWaveDefocusInvertParameterizedARC(550,diffFromOptDist,q1bestAll(subjNumAll==subjNum),q2bestAll(subjNumAll==subjNum),q3bestAll(subjNumAll==subjNum))>380 & ...
                  humanWaveDefocusInvertParameterizedARC(550,diffFromOptDist,q1bestAll(subjNumAll==subjNum),q2bestAll(subjNumAll==subjNum),q3bestAll(subjNumAll==subjNum))<780;
        % CONVERT FROM DEFOCUS AT 550NM TO WAVELENGTH IN FOCUS
        wvInFocusCell{end+1} =  humanWaveDefocusInvertParameterizedARC(550,diffFromOptDist(indGood),q1bestAll(subjNumAll==subjNum),q2bestAll(subjNumAll==subjNum),q3bestAll(subjNumAll==subjNum));
        optDistCnd(end+1,:) = optDistToCheckAll(j);
        rgbLumNormCnd(end+1,:) = conditionsOrderedNorm(i,:);
    end
end

end
