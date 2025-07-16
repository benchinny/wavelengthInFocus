%% LOAD MAIN EXPERIMENT FILES

function [wvInFocusCell, defocusAt550cell, optDistCnd, rgbLumNormCnd] = ARCnlz_lumEffectCheck(subjNum)

bSave = false;
filePath = '/Users/benjaminchin/Documents/ARchromaScraps/meeting_Sept25/';

if subjNum==11
   % blockNums = 2:7;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   blockNums = 11:16;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};   
   subjName = ['S' num2str(subjNum) '-OD'];
elseif subjNum==12
   % blockNums = 2:7;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   % subjName = ['S' num2str(subjNum) '-OD'];

   blockNums = [11:15 18];
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];
elseif subjNum==13
   % blockNums = 3:8;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   blockNums = 12:17;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};   
   subjName = ['S' num2str(subjNum) '-OD'];   
elseif subjNum==14
   blockNums = 9:14;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   % blockNums = 3:8;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};   
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==15
   blockNums = 3:8;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];   
elseif subjNum==17
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==19
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];         
elseif subjNum==20
   blockNums = 3:8;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];         
elseif subjNum==21
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];            
elseif subjNum==22
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];          
elseif subjNum==26
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==27
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==28
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];   
elseif subjNum==30
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
end

meanC = [];
rgb1all = [];
meanv00all = [];
nIndBadTracker = [];
c4all = {};

for i = 1:length(blockNums)
    blockNumTmp = blockNums(i);
    trialNumsTmp = trialNums{i};
    AFCp = ARCloadFileBVAMS(subjNum,blockNumTmp);
    for j = 1:length(trialNumsTmp)
        [ZernikeTable, ~, ~, ~] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumsTmp(j),0);
        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65.
        indBadPupil = table2array(ZernikeTable(:,5))==0;
        PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(~indBadPupil,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0;
        nIndBadTracker(end+1) = sum(indBad);
        c(indBad,4) = mean(c(~indBad,4));
        meanC(end+1,:) = mean(c(1:end,:),1); % TAKE MEAN OF COEFFICIENTS    
        c4all{end+1} = c(:,4);
    end
    rgb1all = [rgb1all; AFCp.rgb100(trialNumsTmp,:)];
    meanv00all = [meanv00all; AFCp.meanv00(trialNumsTmp)./1.2255];
end

% GETTING DEFOCUS AT 550NM
defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
defocusAt550 = humanWaveDefocusARC(550,875,subjNum-10)+meanC(:,4)./defocusCorrectionFactor;

% SORTING CONDITIONS BY COLOR
lumScaleRGB = [4.0888 9.6669 1];

gammaRGB = [2.5 2.7 2.3];

rgbLumNorm = [lumScaleRGB(1).*rgb1all(:,1).^gammaRGB(1) lumScaleRGB(2).*rgb1all(:,2).^gammaRGB(2) lumScaleRGB(3).*rgb1all(:,3).^gammaRGB(3)];

conditionsOrderedNorm = [1.000 0.000 1.000; ...
                         0.625 0.000 0.625];

optDistToCheckAll = [1.5 2.5 3.5];
optDistCnd = [];
rgbLumNormCnd = [];
defocusAt550cell = {};
wvInFocusCell = {};

for j = 1:length(optDistToCheckAll)
    optDistToCheck = optDistToCheckAll(j);
    indDist = abs(meanv00all-optDistToCheck)<0.01;
    for i = 1:size(conditionsOrderedNorm,1)
        ind = abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
              abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
              abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01 & ...
              abs(meanv00all-optDistToCheck)<0.01;
        defocusAt550tmp = defocusAt550(ind);
        diffFromOptDist = defocusAt550tmp-meanv00all(ind);
        indGood = abs(diffFromOptDist)<1 & ...
                  humanWaveDefocusInvertARC(550,diffFromOptDist,subjNum-10)>380 & ...
                  humanWaveDefocusInvertARC(550,diffFromOptDist,subjNum-10)<780;
        wvInFocusCell{end+1} =  humanWaveDefocusInvertARC(550,diffFromOptDist(indGood),subjNum-10);
        defocusAt550cell{end+1} = -defocusAt550tmp(indGood);
        optDistCnd(end+1,:) = optDistToCheckAll(j);
        rgbLumNormCnd(end+1,:) = conditionsOrderedNorm(i,:);
    end
end

end
