%%

%%

filePath = 'G:\My Drive\exp_bvams\code_repo\ARC\';

if strcmp(getenv('username'),'bankslab')
   dataDirectory = filePath;
elseif strcmp(getenv("USER"),'benjaminchin')
   dataDirectory = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Analysis/'; 
elseif strcmp(getenv("USER"),'emily')
   dataDirectory = '/Users/emily/Library/CloudStorage/GoogleDrive-emilyacooper@gmail.com/Shared drives/ARChroma/Analysis/';
end
filenames = {
               [dataDirectory 'S1007V9_AFC_RightACL0_2401181411.mat'] ...
%               [dataDirectory 'S1007V9_AFC_RightACL0_2401181417.mat'] ...
%               [dataDirectory 'S1007V9_AFC_RightACL0_2401181420.mat'] ...  
             };

rgb = [];
meanFocstmOptDst = [];
frqCpd = [];
rspAcu = [];
stimOrientation = [];

for i = 1:length(filenames)
    load(filenames{i});
    rgb = [rgb; AFCp.rgb(1:length(AFCp.rspAcu),:)];
    meanFocstmOptDst = [meanFocstmOptDst; AFCp.meanFocstmOptDst];
    frqCpd = [frqCpd; AFCp.frqCpd((1:length(AFCp.rspAcu)))];
    rspAcu = [rspAcu; AFCp.rspAcu'];
    stimOrientation = [stimOrientation; AFCp.stimOrientation(1:length(AFCp.rspAcu))];
end
