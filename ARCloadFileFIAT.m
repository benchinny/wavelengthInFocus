function [trialData, defocus, strehl, TimeStamp] = ARCloadFileFIAT(subjName,blockNum,trialNum,bExtended,dataPath)

strehl = [];
if ispc
    slash = '\';
else
    slash = '/';
end
dataFolder = [dataPath 'data' slash 'csvFiles' slash];

dirAll = dir(dataFolder);

filenamesAll = {dirAll.name};

filenameTmp = [subjName '-Block_' num2str(blockNum) '-Trial_' num2str(trialNum)];

if bExtended % IF USING EXTENDED ANALYSIS
    ind = find(contains(filenamesAll,filenameTmp) & contains(filenamesAll,'FullPupil'));
else
    ind = find(contains(filenamesAll,filenameTmp) & ~contains(filenamesAll,'FullPupil'));
end

filename = filenamesAll{ind};

trialData = readtable([dataFolder filename]);

if bExtended
   strehl = trialData.Defocus_MaxStrehl;
   defocus = trialData.RMS_LCACorrected_;
else
   defocus = trialData.Defocus_LCACorrected_;
end

TimeStamp = trialData.TimeStamp;

end