function [defocus875,rgbAll,optDistAll,PupilSize,meanCall,q1,q2,q3] = ARCnlzLoadDefocusAbb(subjNum,dataPath)

% function for loading relevant wavefront data from a subject and
% processing it for analysis by other functions (e.g. modeling, weight
% fitting)

% subjNum: subject number. Valid subject numbers: 1, 3, 5, 10, 16, 17, 18,
%          20
% dataPath: local folder where data lives
%
% defocus875: defocus at 875nm
% rgbAll: color conditions
% optDistAll: stimulus distances
% PupilSize: pupil size (for ARC experiment, should always be 4mm)
% meanCall : mean of all Zernike coefficients

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
cAll = []; % MATRIX FOR STORING ALL ZERNIKE COEFFICIENTS
PupilSizeAll = []; % VECTOR FOR STORING PUPIL SIZES PER TRIAL

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
        indBadPupil = table2array(ZernikeTable(:,5))==0; % IDENTIFY BLINKS
        PupilSizeAll(end+1)=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        % ALL ZERNIKE COEFFICIENTS
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        % STACK UP ZERNIKE COEFFICIENTS
        cAll = [cAll; c];
        % FIND FAILURES TO FIT WAVEFRONT IN DEFOCUS VECTOR. NOTE THAT WE
        % ARE DOING THIS FOR A SINGLE TRIAL BECAUSE WE ARE STORING MEAN
        % DEFOCUS SEPARATELY FOR EACH TRIAL. LATER, WE WILL ALSO STORE THE
        % OVERALL AVERAGED WAVEFRONT, SO WE WILL DO THE SAME STEP AGAIN,
        % JUST ACROSS ALL TRIALS.
        indBad = c(:,4)==0;
        meanC = mean(c(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS WITHOUT FAILURES
        % STANDARD CORRECTION TO CONVERT TO EQUIVALENT DEFOCUS
        defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PupilSizeAll(end)/2000)^2);
        defocus875(end+1,:) = meanC(4)./defocusCorrectionFactor;
    end
end

PupilSize = mean(PupilSizeAll); % MEAN PUPIL SIZE (SHOULD BE 4)
% TAKE MEAN OF COEFFICIENTS WITHOUT FAILURES TO FIT WAVEFRONT. NOTE THAT WE
% ARE DOING THIS ACROSS ALL TRIALS, UNLIKE FOR THE VARIABLE 'indBad'
% EARLIER
indBadAll = cAll(:,4)==0; 
meanCall = mean(cAll(~indBadAll,:),1);

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
