function [defocus875rawCell, timeStampRawCell] = ARCmakeFig7Bhelper(subjNum,rgb,optDist,dataPath)

% helper function for grabbing raw accommodative traces from various
% conditions to make Figure 7B

% subjNum: subject number. Valid subject numbers: 1, 3, 5, 10, 16, 17, 18,
%          20
% rgb: RGB values for color condition
% optDist: stimulus distance
% dataPath: local folder where data lives
%
% defocus875rawCell: cell of raw defocus traces
% timeStampRawCell: cell of corresponding time stamps in seconds

% CELL CONTAINING RAW DEFOCUS TRACES FOR CONDITIONS OF INTEREST
defocus875rawCell = {};
timeStampRawCell = {};

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

defocus875rawAll = {}; % CELL FOR STORING ALL RAW DEFOCUS TRACES
optDistAll = []; % STIMULUS DISTANCES
rgbAll = []; % COLOR CONDITIONS
PupilSizeAll = []; % VECTOR FOR STORING PUPIL SIZES PER TRIAL
timeStampSecondsAll = {}; % CELL FOR STORING ALL TIME STAMPS PER TRIAL

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
        % FIND FAILURES TO FIT WAVEFRONT IN DEFOCUS VECTOR. NOTE THAT WE
        % ARE DOING THIS FOR A SINGLE TRIAL BECAUSE WE ARE STORING MEAN
        % DEFOCUS SEPARATELY FOR EACH TRIAL. LATER, WE WILL ALSO STORE THE
        % OVERALL AVERAGED WAVEFRONT, SO WE WILL DO THE SAME STEP AGAIN,
        % JUST ACROSS ALL TRIALS.
        indBad = c(:,4)==0;
        % % REPLACE THEM WITH NANS
        % c(indBad,4) = NaN;
        % OR WITH MEAN ACROSS TRIAL
        c(indBad,4) = mean(c(~indBad,4));
        % STANDARD CORRECTION TO CONVERT TO EQUIVALENT DEFOCUS
        defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PupilSizeAll(end)/2000)^2);
        % CONVERT DEFOCUS ABERRATION COEFFICIENT TO EQUIVALENT DEFOCUS
        defocus875rawAll{end+1} = c(:,4)./defocusCorrectionFactor;
        % CONVERT TIME STAMPS FROM DURATION OBJECT TO SECOND
        timeStampSecondsAll{end+1} = seconds(TimeStamp)-seconds(TimeStamp(1));
    end
end

% INDICES OF TRIALS OF INTEREST
ind = find(abs(rgbAll(:,1)-rgb(1))<0.001 & ...
           abs(rgbAll(:,2)-rgb(2))<0.001 & ...
           abs(rgbAll(:,3)-rgb(3))<0.001 & ...
           abs(optDistAll-optDist)<0.001);

% GRAB AND STORE MEASUREMENTS AND TIMESTAMPS FROM TRIALS OF INTEREST
for i = 1:length(ind)
    defocus875rawCell{end+1} = defocus875rawAll{ind(i)};
    timeStampRawCell{end+1} = timeStampSecondsAll{ind(i)};
end

end