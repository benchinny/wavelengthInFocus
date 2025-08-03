function AFCp = ARCloadFileBVAMS(subjNum,blockNum,dataPath)

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