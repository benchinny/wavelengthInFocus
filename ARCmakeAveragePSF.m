function [meanMeanC,siPSFData] = ARCmakeAveragePSF()

% NOTE SUBJECT NUMBER CONVENTION: SUBTRACT 10 FROM subjNum TO GET ACTUAL
% SUBJECT NUMBER. subjNum VALUES <=10 WERE INTENTIONALLY NOT USED FOR
% ACTUAL PARTICIPANTS. NOTE ALSO THAT PARTICIPANTS WHO DID NOT PASS
% SCREENING OR HAD TO BE EXCLUDED FROM THE ACTUAL ANALYSIS ARE STILL
% INCLUDED IN THIS FUNCTION. 

% subjNum values for participants who passed screening: 11, 13, 15, 20, 26,
% 27, 28, 30. That is, subjects S1, S3, S5, S10, S16, S17, S18, S20. 

%% Initialize and clear
ieInit;

%% Set up display struct and build Ben's stimulus

subjNumAll = [1 3 5 10 16 17 18 20];
meanC = [];

for i = 1:length(subjNumAll)
    subjNum = subjNumAll(i);
    subjNumEncode = subjNum+10;
    
    if subjNum==3
        subjName = 'S13-OD';
        blockNums = 12:17;
        trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
        % blockNums = [2 3];
        % trialNums = [[1:20]' [1:20]']; 
        nTrialTotal = 216;
    elseif subjNum==10
        subjName = 'S20-OD';
        blockNums = 3:8;
        trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
        % blockNums = [2 3];
        % trialNums = [[1:20]' [1:20]'];     
        nTrialTotal = 216;
    elseif subjNum==1
       blockNums = 11:16;
       trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
       subjName = ['S' num2str(subjNum+10) '-OD'];
       nTrialTotal = 216;
    elseif subjNum==5
       blockNums = 3:8;
       trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
       subjName = ['S' num2str(subjNum+10) '-OD'];
       nTrialTotal = 216;   
    elseif subjNum==9
       blockNums = 2:7;
       trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
       subjName = ['S' num2str(subjNum+10) '-OD'];
       nTrialTotal = 216;      
    elseif subjNum==16
       blockNums = 2:7;
       trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
       subjName = ['S' num2str(subjNum+10) '-OD'];
       nTrialTotal = 216;
    elseif subjNum==17
       blockNums = 2:7;
       trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
       subjName = ['S' num2str(subjNum+10) '-OD'];
       nTrialTotal = 216;
    elseif subjNum==18
       blockNums = 2:7;
       trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
       subjName = ['S' num2str(subjNum+10) '-OD'];
       nTrialTotal = 216; 
    elseif subjNum==20
       blockNums = 2:7;
       trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
       subjName = ['S' num2str(subjNum+10) '-OD'];
       nTrialTotal = 216;
    elseif subjNum==21
       blockNums = 2:7;
       trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
       subjName = ['S' num2str(subjNum+9) '-OD'];
       nTrialTotal = 216;   
       subjNumEncode = subjNum+9;
    end
    
    cAll = [];
    optDistAll = [];
    rgbAll = [];
    
    for l = 1:length(blockNums) % LOOP OVER BLOCK
        blockNumInd = l;
        blockNumTmp = blockNums(blockNumInd);
        AFCp = ARCloadFileBVAMS(subjNumEncode,blockNumTmp); % LOAD BVAMS DATA
        rgbAll = [rgbAll; AFCp.rgb100];
        for k = 1:36 % LOOP OVER TRIAL
            trialNumTmp = trialNums(k,blockNumInd);
            
            % LOAD ZERNIKE TABLE AND TIMESTAMPS
            [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumTmp,0);
    
            NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
            c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
            PARAMS = struct;
            PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
            PARAMS.PupilFitSize=mean(table2array(ZernikeTable(:,5))); 
            PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
            c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
            optDistTmp = (AFCp.meanv00(k)./1.2255).*ones([size(c,1) 1]);
            optDistAll = [optDistAll; optDistTmp];
            cAll = [cAll; c];
        end
    end
    
    indBad = cAll(:,4)==0;
    meanC(i,:) = mean(cAll(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS
end

meanMeanC = 1.*mean(meanC,1);

wave = 380:4:780;
zCoeffs = [0 meanMeanC(1:end-1)];
wvfP = wvfCreate('calc wavelengths', wave, ...
    'measured wavelength', 552, ...
    'zcoeffs', zCoeffs, 'measured pupil', PARAMS.PupilSize, ...
    'name', sprintf('human-%d', PARAMS.PupilSize),'spatial samples',260);
wvfP.calcpupilMM = PARAMS.PupilSize;
wvfP.refSizeOfFieldMM = 12;

wvfP = wvfSet(wvfP, 'zcoeff', 0, 'defocus');
        
% Convert to siData format as well as wavefront object
[siPSFData, wvfP] = wvf2SiPsfARC(wvfP,'showBar',false,'nPSFSamples',260,'umPerSample',1.1512); 
end
