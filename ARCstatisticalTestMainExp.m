%% JUST SETTING UP PATHS TO PRE-SAVED DATA

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end

load([dataPath 'data' slash 'PresavedFigureData' slash 'allExp1DataRGB.mat']);

%% CREATE COLUMNS FOR TABLE FOR LINEAR MODEL

% CODING EACH PARTICIPANT NUMBER AS A LETTER SO MATLAB DOESN'T INTERPRET
% PARTICIPANT NUMBERS AS CONTINUOUS NUMERIC DATA
subjStr = 'A B C    D     EFG H';

rgbLumNormCndUnq = unique(rgbLumNormCndAll,'rows'); % UNIQUE CONDITIONS
indGreenless = [1 3 7 6 5]; % INDICES FOR ORDERING BLUEST TO REDDEST
indGreen = [2 4 10 9 8]; % SAME AS ABOVE BUT FOR 'SOME GREEN' CONDITIONS

rbRatioGreenless = []; % RED BLUE RATIO
optDistCndGreenless = []; % STIMULUS DISTANCE
subjNumTagGreenless = []; % SUBJECT NUMBER TAG
wvInFocusMeanGreenless = []; % WAVELENGTH IN FOCUS

% THIS LOOP IS JUST SORTING DATA INTO COLUMNS OF TABLE--SEE 'dataTable'
% VARIABLES IN NEXT SECTION
for i = 1:length(indGreenless) 
    % FOR EACH COLOR, GET THE MEAN ACCOMMODATIVE RESPONSE FOR EACH
    % PARTICIPANT AND OPTICAL DISTANCE
    rgbTmp = rgbLumNormCndUnq(indGreenless(i),:); 
    indRgbTmp = find(abs(rgbLumNormCndAll(:,1)-rgbTmp(1))<0.001 & ... 
                     abs(rgbLumNormCndAll(:,2)-rgbTmp(2))<0.001 & ... 
                     abs(rgbLumNormCndAll(:,3)-rgbTmp(3))<0.001);
    for j = 1:length(indRgbTmp) % ADD TO VECTOR
        rbRatioGreenless(end+1) = rgbLumNormCndAll(indRgbTmp(j),1)/rgbLumNormCndAll(indRgbTmp(j),3);
        optDistCndGreenless(end+1) = optDistCndAll(indRgbTmp(j));
        subjNumTagGreenless = [subjNumTagGreenless; subjStr(subjNumTag(indRgbTmp(j)))];
        wvInFocusMeanGreenless(end+1) = mean(wvInFocusCellAll{indRgbTmp(j)});
    end
end

% SAME AS ABOVE LOOP, JUST FOR STIMULI WITH GREEN PRIMARY ON
rbRatioGreen = [];
optDistCndGreen = [];
subjNumTagGreen = [];
wvInFocusMeanGreen = [];

for i = 1:length(indGreen) 
    rgbTmp = rgbLumNormCndUnq(indGreen(i),:); 
    indRgbTmp = find(abs(rgbLumNormCndAll(:,1)-rgbTmp(1))<0.001 & ... 
                     abs(rgbLumNormCndAll(:,2)-rgbTmp(2))<0.001 & ... 
                     abs(rgbLumNormCndAll(:,3)-rgbTmp(3))<0.001);
    for j = 1:length(indRgbTmp)
        rbRatioGreen(end+1) = rgbLumNormCndAll(indRgbTmp(j),1)/rgbLumNormCndAll(indRgbTmp(j),3);
        optDistCndGreen(end+1) = optDistCndAll(indRgbTmp(j));
        subjNumTagGreen = [subjNumTagGreen; subjStr(subjNumTag(indRgbTmp(j)))];
        wvInFocusMeanGreen(end+1) = mean(wvInFocusCellAll{indRgbTmp(j)});
    end
end

%% CREATE TABLES

dataTableGreen = array2table([optDistCndGreen' wvInFocusMeanGreen'], ...
                             'VariableNames',{'OpticalDistance' 'MeanD'});

dataTableGreen.("Subject") = subjNumTagGreen;
dataTableGreen.("ColorRatio") = rbRatioGreen';

dataTableGreenless = array2table([optDistCndGreenless' wvInFocusMeanGreenless'], ...
                             'VariableNames',{'OpticalDistance' 'MeanD'});

dataTableGreenless.("Subject") = subjNumTagGreenless;
dataTableGreenless.("ColorRatio") = rbRatioGreenless';

%% RUN GLME

glmeGreenless = fitglme(dataTableGreenless,'MeanD ~ ColorRatio + OpticalDistance + ColorRatio:OpticalDistance + (1|Subject)', ...
    'Distribution','Normal','Link','identity')

glmeGreen = fitglme(dataTableGreen,'MeanD ~ ColorRatio + OpticalDistance + ColorRatio:OpticalDistance + (1|Subject)', ...
    'Distribution','Normal','Link','identity')

