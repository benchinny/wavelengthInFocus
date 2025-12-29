function AFCp = ARCcopyFileBVAMS(subjNum,blockNum,dataPath)

% This function copies BVAMS data to a new folder

% INITIALIZE SOMETHING TO RETURN EVEN IF NOTHING IS FOUND
AFCp = []; 

dataFolder = fullfile(dataPath,'data','ARC');
dataFolderNew = fullfile(dataPath,'data','data','stimulusData');

% LOAD 'MASTER FILE' SHOWING ALL FILE NAMES FOR EACH BLOCK AND CONDITION
load(fullfile(dataFolder,'AFCflsR'));

% FIND THE FILE NAME FOR THE BLOCK AND CONDITION IN QUESTION
filenameTmp = AFCfls{subjNum,blockNum};

% GRAB THE PART OF THE FILE NAME THAT DOESN'T INCLUDE THE PATH ON THE BVAMS
% LOCAL MACHINE
filenameTmp = filenameTmp((length('H:\Shared drives\CIVO_BVAMS\data\ARC\')+1):end);

copyfile(fullfile(dataFolder,[filenameTmp '.mat']),fullfile(dataFolderNew,[filenameTmp '.mat']));

end