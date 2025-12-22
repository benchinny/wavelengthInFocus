%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

% LIST OF ALL SUBJECTS
subjNumListAll = [1 3 5 10 16 17 18 20];

for i = 1:length(subjNumListAll)
    subjNum = subjNumListAll(i);
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
    
    trialNumAll = 1:36;
    % LOAD DATA TO FIT
    for k = 1:length(blockNumAll) % LOOP OVER BLOCKS
        for l = 1:length(trialNumAll) % FOR EACH TRIAL
           [ZernikeTable, ~, TimeStamp] = ARCcopyFileFIAT(subjName,blockNumAll(k),trialNumAll(l),dataPath);
        end
    end
end

%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
dataFolder = fullfile(dataPath,'data','ARC');
dataFolderNew = fullfile(dataPath,'data','data','psychophysicalData');

% LIST OF ALL SUBJECTS
subjNumListAll = [1 3 5 10 16 17 18 20];

for i = 1:length(subjNumListAll)
    subjNum = subjNumListAll(i);
    if subjNum==3
        filenames = {
                      % % S3 PURPLE (ACTUAL)
                      'S1013V19_AFC_RightACL0_2409201559.mat' ...
                      'S1013V19_AFC_RightACL0_2409201548.mat' ...
                      'S1013V19_AFC_RightACL0_2409201542.mat' ...
                      'S1013V19_AFC_RightACL0_2409201529.mat' ...   
                      };
    elseif subjNum==1
        filenames = {
                      % % S1 PURPLE (ACTUAL)
                      'S1011V18_AFC_RightACL0_2411050943.mat' ...
                      'S1011V18_AFC_RightACL0_2411050928.mat' ...
                      };
    elseif subjNum==10
        filenames = {
                      % S10 PURPLE
                      'S1020V10_AFC_RightACL0_2409111552.mat' ...
                      'S1020V10_AFC_RightACL0_2409111600.mat' ...
                      'S1020V10_AFC_RightACL0_2409111607.mat' ...
                      'S1020V10_AFC_RightACL0_2409111615.mat' ...
                      };
    elseif subjNum==5
        filenames = {
                      % S10 PURPLE
                      'S1015V10_AFC_RightACL0_2410031638.mat' ...
                      'S1015V10_AFC_RightACL0_2410031631.mat' ...
                      'S1015V10_AFC_RightACL0_2410031615.mat' ...
                      'S1015V10_AFC_RightACL0_2410031608.mat' ...
                      };    
    elseif subjNum==16
        filenames = {
                      % S16 PURPLE
                      'S1026V8_AFC_RightACL0_2409231008.mat' ...
                      'S1026V8_AFC_RightACL0_2409231000.mat' ...
                      'S1026V8_AFC_RightACL0_2409230953.mat' ...
                      'S1026V8_AFC_RightACL0_2409230946.mat' ...    
                      };    
    elseif subjNum==17
        filenames = {
                      % S17 PURPLE
                      'S1027V9_AFC_RightACL0_2410081135.mat' ...
                      'S1027V9_AFC_RightACL0_2410081141.mat' ...
                      'S1027V9_AFC_RightACL0_2410081208.mat' ...
                      'S1027V9_AFC_RightACL0_2410081214.mat' ...    
                      }; 
    elseif subjNum==18
        filenames = {
                      % S17 PURPLE
                      'S1028V9_AFC_RightACL0_2411051411.mat' ...
                      'S1028V9_AFC_RightACL0_2411051420.mat' ...
                      'S1028V9_AFC_RightACL0_2411151223.mat' ...
                      'S1028V9_AFC_RightACL0_2411151217.mat' ...    
                      }; 
    elseif subjNum==20
        filenames = {
                      % S17 PURPLE
                      'S1030V8_AFC_RightACL0_2411151125.mat' ...
                      'S1030V8_AFC_RightACL0_2411151134.mat' ...
                      'S1030V8_AFC_RightACL0_2411151143.mat' ...
                      'S1030V8_AFC_RightACL0_2411151152.mat' ...    
                      };         
    end
    for j = 1:length(filenames)
        copyfile(fullfile(dataFolder,[filenames{j}]),fullfile(dataFolderNew,[filenames{j}]));
    end
end