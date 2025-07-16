function [trialDataAll, defocus, strehl, TimeStamp, trialTagIndex] = ARCloadFileFIATallInstances(subjName,blockNum,bExtended)

strehl = [];

if strcmp(getenv("USER"),'ben')
   dataFolder= '/home/ben/Documents/ARchroma/FIAT/csvFiles/';
else
   dataFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/csvFiles/';
end

dirAll = dir(dataFolder);

filenamesAll = {dirAll.name};

filenameTmp = [subjName '-Block_' num2str(blockNum)];

if bExtended % IF USING EXTENDED ANALYSIS
    ind = find(contains(filenamesAll,filenameTmp) & contains(filenamesAll,'FullPupil'));
else
    ind = find(contains(filenamesAll,filenameTmp) & ~contains(filenamesAll,'FullPupil'));
end

trialTagIndex = [];
trialDataAll = table;
for subInd = 1:length(ind)
   filename = filenamesAll{ind(subInd)};
   trialData = readtable([dataFolder filename]);
   trialTagIndex = [trialTagIndex; subInd.*ones([size(trialData,1) 1])];
   if subInd==1
      trialDataAll = trialData;
   else
      trialDataAll = vertcat(trialDataAll,trialData);
   end
end

if bExtended
   strehl = trialData.Defocus_MaxStrehl;
   defocus = trialData.RMS_LCACorrected_;
else
   defocus = trialData.Defocus_LCACorrected_;
end

TimeStamp = trialData.TimeStamp;

end