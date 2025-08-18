function AFCp = ARCloadFileBVAMS(subjNum,blockNum,dataPath)

% THIS FUNCTION LOADS DATA FROM THE BVAMS (.mat FILES)--EVERYTHING THAT
% HAPPENS HERE IS SPECIFIC TO THE NAMING AND ORGANIZATION CONVENTIONS OF
% THE DATA

AFCp = [];
if ispc
    slash = '\';
else
    slash = '/';
end
dataFolder = [dataPath 'data' slash 'ARC' slash];

load([dataFolder 'AFCflsR']);

filenameTmp = AFCfls{subjNum,blockNum};

filenameTmp = filenameTmp((length('H:\Shared drives\CIVO_BVAMS\data\ARC\')+1):end);

load([dataFolder filenameTmp]);

end