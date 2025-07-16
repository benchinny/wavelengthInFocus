%%

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/meeting_Sept25/allExp1DataRGB.mat');

%%

subjStr = 'A B C    D     EFG H';

rgbLumNormCndUnq = unique(rgbLumNormCndAll,'rows');
indGreenless = [1 3 7 6 5];
indGreen = [2 4 10 9 8];

defocusAt550meanGreenless = [];
rbRatioGreenless = [];
optDistCndGreenless = [];
subjNumTagGreenless = [];
wvInFocusMeanGreenless = [];

for i = 1:length(indGreenless) 
    rgbTmp = rgbLumNormCndUnq(indGreenless(i),:); 
    indRgbTmp = find(abs(rgbLumNormCndAll(:,1)-rgbTmp(1))<0.001 & ... 
                     abs(rgbLumNormCndAll(:,2)-rgbTmp(2))<0.001 & ... 
                     abs(rgbLumNormCndAll(:,3)-rgbTmp(3))<0.001);
    for j = 1:length(indRgbTmp)
        defocusAt550meanGreenless(end+1) = mean(defocusAt875cellAll{indRgbTmp(j)});
        rbRatioGreenless(end+1) = rgbLumNormCndAll(indRgbTmp(j),1)/rgbLumNormCndAll(indRgbTmp(j),3);
        optDistCndGreenless(end+1) = optDistCndAll(indRgbTmp(j));
        subjNumTagGreenless = [subjNumTagGreenless; subjStr(subjNumTag(indRgbTmp(j)))];
        wvInFocusMeanGreenless(end+1) = mean(wvInFocusCellAll{indRgbTmp(j)});
    end
end

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

%%

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

%%

glmeGreenless = fitglme(dataTableGreenless,'MeanD ~ ColorRatio + OpticalDistance + ColorRatio:OpticalDistance + (1|Subject)', ...
    'Distribution','Normal','Link','identity');

glmeGreen = fitglme(dataTableGreen,'MeanD ~ ColorRatio + OpticalDistance + ColorRatio:OpticalDistance + (1|Subject)', ...
    'Distribution','Normal','Link','identity');

