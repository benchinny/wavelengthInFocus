function AFCp = ARCloadFileBVAMS(subjNum,blockNum)

AFCp = [];

if strcmp(getenv("USER"),'ben')
   dataFolder= '/home/ben/Documents/ARchroma/FIAT/ARC/';
elseif strcmp(getenv("USER"),'benchin')
   dataFolder = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/ARC/';
else
   dataFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/ARC/';
end

load([dataFolder 'AFCflsR']);

filenameTmp = AFCfls{subjNum,blockNum};

filenameTmp = filenameTmp((length('H:\Shared drives\CIVO_BVAMS\data\ARC\')+1):end);

load([dataFolder filenameTmp]);

end