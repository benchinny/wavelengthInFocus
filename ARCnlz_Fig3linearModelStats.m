%% JUST SETTING UP PATHS TO PRE-SAVED DATA

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

load(fullfile(dataPath,'data','PresavedFigureData','allExp1DataRGB.mat'));

%% CREATE COLUMNS FOR TABLE FOR LINEAR MODEL

subjNumTagOptions = 'ABCDEFGH';
% INDICES OF PRESAVED DATA THAT CORRESPOND TO CONDITIONS WITH NO GREEN IN
% THEM
indGreenless = 1:5;
% INDICES OF PRESAVED DATA THAT CORRESPOND TO CONDITIONS WITH SOME GREEN IN
% THEM
indGreen = 6:10;

% FIRST MAKE APPROPRIATE COLUMN VECTORS FOR CONDITIONS WITH SOME GREEN IN
% THEM

optDistCndGreen = []; % STIMULUS DISTANCES
wvInFocusMeanGreen = []; % WAVELENGTHS IN FOCUS
subjNumTagGreen = []; % ALPHABETICAL TAG FOR EACH SUBJECT
rbRatioGreen = []; % RED-BLUE RATIO

for i = 1:length(subjNumTagOptions) % LOOP OVER SUBJECTS
    for j = 1:size(optDistCndAll,1) % LOOP OVER DISTANCES
        wvInFocusMeanGreen = [wvInFocusMeanGreen; squeeze(wvMeanAll(indGreen,j,i))];
        % GRAB COLOR CONDITIONS FIRST
        conditionsOrderedNormTmp = squeeze(conditionsOrderedNormAll(indGreen,:,i));
        rbRatioGreen = [rbRatioGreen; conditionsOrderedNormTmp(:,1)./conditionsOrderedNormTmp(:,3)];
        optDistCndGreen = [optDistCndGreen; optDistCndAll(j,i).*ones([length(indGreen) 1])];
        subjNumTagGreen = [subjNumTagGreen; repmat(subjNumTagOptions(i),[length(indGreen) 1])];
    end
end

% MAKE APPROPRIATE COLUMN VECTORS FOR CONDITIONS WITH NO GREEN IN THEM

optDistCndGreenless = []; % STIMULUS DISTANCES
wvInFocusMeanGreenless = []; % WAVELENGTHS IN FOCUS
subjNumTagGreenless = []; % ALPHABETICAL TAG FOR EACH SUBJECT
rbRatioGreenless = []; % RED-BLUE RATIO

for i = 1:length(subjNumTagOptions) % LOOP OVER SUBJECTS
    for j = 1:size(optDistCndAll,1) % LOOP OVER DISTANCES
        wvInFocusMeanGreenless = [wvInFocusMeanGreenless; squeeze(wvMeanAll(indGreenless,j,i))];
        % GRAB COLOR CONDITIONS FIRST
        conditionsOrderedNormTmp = squeeze(conditionsOrderedNormAll(indGreenless,:,i));
        rbRatioGreenless = [rbRatioGreenless; conditionsOrderedNormTmp(:,1)./conditionsOrderedNormTmp(:,3)];
        optDistCndGreenless = [optDistCndGreenless; optDistCndAll(j,i).*ones([length(indGreenless) 1])];
        subjNumTagGreenless = [subjNumTagGreenless; repmat(subjNumTagOptions(i),[length(indGreenless) 1])];
    end
end

%% CREATE TABLES

dataTableGreen = array2table([-optDistCndGreen wvInFocusMeanGreen], ...
                             'VariableNames',{'OpticalDistance' 'MeanD'});

dataTableGreen.("Subject") = subjNumTagGreen;
dataTableGreen.("ColorRatio") = rbRatioGreen;

dataTableGreenless = array2table([-optDistCndGreenless wvInFocusMeanGreenless], ...
                             'VariableNames',{'OpticalDistance' 'MeanD'});

dataTableGreenless.("Subject") = subjNumTagGreenless;
dataTableGreenless.("ColorRatio") = rbRatioGreenless;

%% RUN GLME

glmeGreenless = fitglme(dataTableGreenless,'MeanD ~ ColorRatio + OpticalDistance + ColorRatio:OpticalDistance + (1|Subject)', ...
    'Distribution','Normal','Link','identity')

glmeGreen = fitglme(dataTableGreen,'MeanD ~ ColorRatio + OpticalDistance + ColorRatio:OpticalDistance + (1|Subject)', ...
    'Distribution','Normal','Link','identity')

