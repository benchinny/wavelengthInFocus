%% LOAD MAIN EXPERIMENT FILES

function [wvInFocusCell, defocusAt550cell, defocusAt875cell, optDistCnd, rgbLumNormCnd] = ARCnlz_mainExpSortColorVisStrehl(subjNum)

% NOTE SUBJECT NUMBER CONVENTION: SUBTRACT 10 FROM subjNum TO GET ACTUAL
% SUBJECT NUMBER. subjNum VALUES <=10 WERE INTENTIONALLY NOT USED FOR
% ACTUAL PARTICIPANTS. NOTE ALSO THAT PARTICIPANTS WHO DID NOT PASS
% SCREENING OR HAD TO BE EXCLUDED FROM THE ACTUAL ANALYSIS ARE STILL
% INCLUDED IN THIS FUNCTION. 

% subjNum values for participants who passed screening: 11, 13, 15, 20, 26,
% 27, 28, 30. That is, subjects S1, S3, S5, S10, S16, S17, S18, S20. 

% This function grabs all raw wavefront data from Experiment 1 for a
% particular subject and returns mean defocus values (both at 550nm and 875nm),
% wavelength-in-focus values, and the optical distances and color
% conditions corresponding to those values. The 'c' variable contains the
% raw traces of every coefficient on each of the 65 Zernike terms. 

bSave = false;
bStrehl = true;
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
rgb1all = []; % RGB VALUES
meanv00all = []; % OPTOTUNE VALUE (TO GET STIMULUS DISTANCE, DIVIDE BY 1.2255)
nIndBadTracker = [];
c4all = {};

% GRABBING ZERNIKE VALUES FROM DATA REPOSITORY
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
        indBad = c(:,4)==0; % VALUE OF EXACTLY 0 MEANS BLINK
        nIndBadTracker(end+1) = sum(indBad);
        c(indBad,4) = mean(c(~indBad,4));
        meanC(end+1,:) = mean(c(1:end,:),1); % TAKE MEAN OF COEFFICIENTS    
        c4all{end+1} = c(:,4); % COEFFICIENTS ON DEFOCUS TERM
    end
    rgb1all = [rgb1all; AFCp.rgb100(trialNumsTmp,:)];
    meanv00all = [meanv00all; AFCp.meanv00(trialNumsTmp)./1.2255];
end

% GETTING DEFOCUS AT 550NM
defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
defocusAt550 = humanWaveDefocusARC(550,875,subjNum-10)+meanC(:,4)./defocusCorrectionFactor;
defocusAt875 = meanC(:,4)./defocusCorrectionFactor;

% SORTING CONDITIONS BY COLOR
lumScaleRGB = [4.0888 9.6669 1];

gammaRGB = [2.5 2.7 2.3];

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

optDistToCheckAll = [1.5 2.5 3.5];
optDistCnd = [];
rgbLumNormCnd = [];
defocusAt550cell = {};
defocusAt875cell = {};
wvInFocusCell = {};

% % % making 2D CSF function
f1Dhalf = 0.8105:0.8105:(0.8105*159);
f1D = [-fliplr(f1Dhalf) 0 f1Dhalf(1:end-1)];
[fx, fy] = meshgrid(f1D,f1D);
df = sqrt(fx.^2 + fy.^2); % compute distance from origin
CSF2d = 0.04992*(1+5.9375*df).*exp(-0.114*df.^1.1);
% inverse Fourier transform of 2D CSF
N = ifftshift(ifft2(fftshift(CSF2d)));

wvfPdiffLim = wvfCreate('calc wavelengths', 875, ...
                        'measured wavelength', 875, ...
                        'zcoeffs', zeros([1 65]), 'measured pupil', PARAMS.PupilSize, ...
                        'name', sprintf('human-%d', PARAMS.PupilSize),'spatial samples',318);
wvfPdiffLim.calcpupilMM = PARAMS.PupilSize;
wvfPdiffLim.refSizeOfFieldMM = 12;
wvfPdiffLim = wvfSet(wvfPdiffLim, 'zcoeff', 0, 'defocus');
[siPSFDataDiffLim, wvfPdiffLim] = wvf2SiPsfARC(wvfPdiffLim,'showBar',false,'nPSFSamples',318,'umPerSample',1.1512); 
psfDiffLimWeighted = siPSFDataDiffLim.psf.*N;
vsxDenom = sum(psfDiffLimWeighted(:));

for j = 1:length(optDistToCheckAll)
    optDistToCheck = optDistToCheckAll(j);
    indDist = abs(meanv00all-optDistToCheck)<0.01;
    for i = 1:size(conditionsOrderedNorm,1)
        ind = abs(rgbLumNorm(:,1)-conditionsOrderedNorm(i,1))<0.01 & ...
              abs(rgbLumNorm(:,2)-conditionsOrderedNorm(i,2))<0.01 & ...
              abs(rgbLumNorm(:,3)-conditionsOrderedNorm(i,3))<0.01 & ...
              abs(meanv00all-optDistToCheck)<0.01;
        if bStrehl
            meanCtmp = mean(meanC(ind,:),1);
            zCoeffs = [0 meanCtmp(1:end-1)];
            defocusToTest = -0.5:0.05:0.5;
            peakPSF = [];
            for k = 1:length(defocusToTest)
                wvfP = wvfCreate('calc wavelengths', 875, ...
                    'measured wavelength', 875, ...
                    'zcoeffs', zCoeffs, 'measured pupil', PARAMS.PupilSize, ...
                    'name', sprintf('human-%d', PARAMS.PupilSize),'spatial samples',318);
                wvfP.calcpupilMM = PARAMS.PupilSize;
                wvfP.refSizeOfFieldMM = 12;
                wvfP = wvfSet(wvfP, 'zcoeff', defocusToTest(k)*defocusCorrectionFactor, 'defocus');
                [siPSFData, wvfP] = wvf2SiPsfARC(wvfP,'showBar',false,'nPSFSamples',318,'umPerSample',1.1512); 
                psfWeighted = siPSFData.psf.*N;
                vsxNum = sum(psfWeighted(:));
                peakPSF(k) = vsxNum./vsxDenom;
                display(['Iterate through defocus: ' num2str(defocusToTest(k))]);
            end
            [~,indMaxStrehl] = max(peakPSF);
            defocusOffset = defocusToTest(indMaxStrehl);
            defocusAt550tmp = defocusAt550(ind)-defocusOffset;
            defocusAt875tmp = defocusAt875(ind)-defocusOffset;
            diffFromOptDist = defocusAt550tmp-meanv00all(ind);
            % EXCLUDE DATA FOR WHICH PARTICIPANT WAS ACCOMMODATING OUTSIDE OF
            % VISIBLE RANGE
            indGood = abs(diffFromOptDist)<1 & ...
                      humanWaveDefocusInvertARC(550,diffFromOptDist,subjNum-10)>380 & ...
                      humanWaveDefocusInvertARC(550,diffFromOptDist,subjNum-10)<780;
            wvInFocusCell{end+1} =  humanWaveDefocusInvertARC(550,diffFromOptDist(indGood),subjNum-10);
            defocusAt550cell{end+1} = -defocusAt550tmp(indGood);
            defocusAt875cell{end+1} = -defocusAt875tmp(indGood);
            optDistCnd(end+1,:) = optDistToCheckAll(j);
            rgbLumNormCnd(end+1,:) = conditionsOrderedNorm(i,:);
        end
    end
end

end
