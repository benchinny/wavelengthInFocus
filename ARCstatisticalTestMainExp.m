%%

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end

load([dataPath 'data' slash 'PresavedFigureData' slash 'allExp1DataRGB.mat']);

%%

% CODING EACH PARTICIPANT NUMBER AS A LETTER SO MATLAB DOESN'T INTERPRET
% PARTICIPANT NUMBERS AS CONTINUOUS NUMERIC DATA
subjStr = 'A B C    D     EFG H';

rgbLumNormCndUnq = unique(rgbLumNormCndAll,'rows'); % UNIQUE CONDITIONS
indGreenless = [1 3 7 6 5];
indGreen = [2 4 10 9 8];

% ORIGINALLY THIS VARIABLE WAS DEFOCUS AT 550, BUT I CHANGED IT TO DEFOCUS
% AT 875 WITHOUT CHANGING THIS VARIABLE NAME. I SHOULD FIX THIS. 
defocusAt550meanGreenless = [];
rbRatioGreenless = [];
optDistCndGreenless = [];
subjNumTagGreenless = [];
wvInFocusMeanGreenless = [];

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
        defocusAt550meanGreenless(end+1) = mean(defocusAt875cellAll{indRgbTmp(j)});
        rbRatioGreenless(end+1) = rgbLumNormCndAll(indRgbTmp(j),1)/rgbLumNormCndAll(indRgbTmp(j),3);
        optDistCndGreenless(end+1) = optDistCndAll(indRgbTmp(j));
        subjNumTagGreenless = [subjNumTagGreenless; subjStr(subjNumTag(indRgbTmp(j)))];
        wvInFocusMeanGreenless(end+1) = mean(wvInFocusCellAll{indRgbTmp(j)});
    end
end

% SAME AS ABOVE LOOP, JUST FOR STIMULI WITH GREEN PRIMARY ON
defocusAt550meanGreen = [];
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
        defocusAt550meanGreen(end+1) = mean(defocusAt875cellAll{indRgbTmp(j)});
        rbRatioGreen(end+1) = rgbLumNormCndAll(indRgbTmp(j),1)/rgbLumNormCndAll(indRgbTmp(j),3);
        optDistCndGreen(end+1) = optDistCndAll(indRgbTmp(j));
        subjNumTagGreen = [subjNumTagGreen; subjStr(subjNumTag(indRgbTmp(j)))];
        wvInFocusMeanGreen(end+1) = mean(wvInFocusCellAll{indRgbTmp(j)});
    end
end

% rbRatioGreen = ordinal(rbRatioGreen);
% rbRatioGreenless = ordinal(rbRatioGreenless);

%% CREATE TABLES

% defocusAt550meanGreen = defocusAt550meanGreen(randperm(120));

% defocusAt550meanGreenless = defocusAt550meanGreenless(randperm(120));

% rbRatioGreen = rbRatioGreen(randperm(120));
% 
% rbRatioGreenless = rbRatioGreenless(randperm(120));

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
    'Distribution','Normal','Link','identity');

glmeGreen = fitglme(dataTableGreen,'MeanD ~ ColorRatio + OpticalDistance + ColorRatio:OpticalDistance + (1|Subject)', ...
    'Distribution','Normal','Link','identity');

