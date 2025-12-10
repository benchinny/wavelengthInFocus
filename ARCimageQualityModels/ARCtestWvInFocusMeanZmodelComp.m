%% FIGURE 4C

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];

% LOAD PRE-SAVED DATA AND MODEL FITS
load([foldername 'wvMeanAndPredLM.mat']);
aicSpatFilterLM = aicAll;
load([foldername 'wvMeanAndPredLminusM.mat']);
aicSpatFilterLminusM = aicAll;
load([foldername 'wvMeanAndPredDonutx2.mat']);
aicSpatFilterLMS = aicAll;

figure; 
hold on;
boxplot([(aicSpatFilterLM-aicSpatFilterLMS)' ...
          (aicSpatFilterLM-aicSpatFilterLminusM)']);
plot(1,(aicSpatFilterLM-aicSpatFilterLMS)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
plot(2,(aicSpatFilterLM-aicSpatFilterLminusM)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
ylim([-25 120]);
set(gca,'YTick',[-23.0259 0 23.0259 46.0517 69.0776 92.1034 115.1293]);
set(gca,'YTickLabel',{'10^-10' '1' '10^10' '10^20' '10^30' '10^40' '10^50'});
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Blue-yellow vs Luminance' 'Red-green vs Luminance'});

%% FIGURE 4C RESPONSE TO REVIEWER WITH STREHL

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];

% LOAD PRE-SAVED DATA AND MODEL FITS
load([foldername 'wvMeanAndPredStrehlLM.mat']);
aicStrehlFilterLM = aicAll;
load([foldername 'wvMeanAndPredStrehlLminusM.mat']);
aicStrehlLminusM = aicAll;
load([foldername 'wvMeanAndPredStrehlLMS.mat']);
aicStrehlLMS = aicAll;

figure; 
hold on;
boxplot([(aicStrehlFilterLM-aicStrehlLMS)' ...
          (aicStrehlFilterLM-aicStrehlLminusM)']);
plot(1,(aicStrehlFilterLM-aicStrehlLMS)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
plot(2,(aicStrehlFilterLM-aicStrehlLminusM)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
ylim([-25 120]);
set(gca,'YTick',[-23.0259 0 23.0259 46.0517 69.0776 92.1034 115.1293]);
set(gca,'YTickLabel',{'10^-10' '1' '10^10' '10^20' '10^30' '10^40' '10^50'});
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Blue-yellow vs Luminance' 'Red-green vs Luminance'});

%% SUPPLEMENTAL FIGURE WITH FINCH MODEL

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];

% LOAD PRE-SAVED DATA AND MODEL FITS
load([foldername 'wvMeanAndPredDeltaPass2.mat']);
aicSpatFilterLM = aicAll;
load([foldername 'wvMeanAndPredLminusM.mat']);
aicSpatFilterLminusM = aicAll;
load([foldername 'wvMeanAndPredDonutx2.mat']);
aicSpatFilterLMS = aicAll;

figure; 
hold on;
boxplot([(aicSpatFilterLM-aicSpatFilterLMS)' ...
          (aicSpatFilterLM-aicSpatFilterLminusM)']);
plot(1,(aicSpatFilterLM-aicSpatFilterLMS)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
plot(2,(aicSpatFilterLM-aicSpatFilterLminusM)','k.','MarkerSize',10,'MarkerFaceColor',[0 0 0]);
ylim([-25 120]);
set(gca,'YTick',[-23.0259 0 23.0259 46.0517 69.0776 92.1034 115.1293]);
set(gca,'YTickLabel',{'10^-10' '1' '10^10' '10^20' '10^30' '10^40' '10^50'});
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Blue-yellow vs Luminance' 'Red-green vs Luminance'});
