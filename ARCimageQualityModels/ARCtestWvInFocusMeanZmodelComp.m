%% SPECIFY PATH TO DATA FOLDER

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

%% LOADING AIC VALUES AND USING THEM TO MAKE PLOTS

if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'AICmodelComparisons' slash];

load([foldername 'aicDeltaPass2.mat'])
load([foldername 'aicDeltaPass4.mat'])
load([foldername 'aicSpatFilterLM.mat'])
load([foldername 'aicSpatFilterLminusM.mat'])
load([foldername 'aicSpatFilterLMS.mat'])
load([foldername 'aicStrehlLMS.mat'])
load([foldername 'aicStrehlLminusM.mat'])
load([foldername 'aicStrehlLM.mat'])

barWidth = 0.1;

figure;
set(gcf,'Position',[472 453 1054 420]);
hold on;
bar(1:10:71,aicSpatFilterLMS,barWidth,'FaceColor',[0 0 1]);
bar(2:10:72,aicSpatFilterLminusM,barWidth,'FaceColor',[1 0.5 0]);
bar(3:10:73,aicSpatFilterLM,barWidth,'FaceColor','w');
bar(4:10:74,aicStrehlLMS,barWidth);
bar(5:10:75,aicDeltaPass2,barWidth,'FaceColor',0.7.*[1 1 1]);
bar(6:10:76,aicDeltaPass4,barWidth,'FaceColor',0.4.*[1 1 1]);
bar(7:10:77,aicStrehlLminusM,barWidth,'FaceColor',1.*[0.7 0.45 0]);
bar(8:10:78,aicStrehlLM,barWidth,'FaceColor',1.*[1 1 1]);
plot(-1.*[1 1],[-100 50],'k-','LineWidth',1);
plot(9.*[1 1],[-100 50],'k-','LineWidth',1);
plot(19.*[1 1],[-100 50],'k-','LineWidth',1);
plot(29.*[1 1],[-100 50],'k-','LineWidth',1);
plot(39.*[1 1],[-100 50],'k-','LineWidth',1);
plot(49.*[1 1],[-100 50],'k-','LineWidth',1);
plot(59.*[1 1],[-100 50],'k-','LineWidth',1);
plot(69.*[1 1],[-100 50],'k-','LineWidth',1);
plot(79.*[1 1],[-100 50],'k-','LineWidth',1);
xlim([-2 80]);
set(gca,'XTick',3.5:9:66.5);
set(gca,'FontSize',15);
set(gca,'XTickLabel',{'S1' 'S2' 'S3' 'S4' 'S5' 'S6' 'S7' 'S8'})
ylabel('AIC');
legend('L+M-S','L-M','L+M','L+M-S Strehl','L+M 2cpd','L+M 4cpd','L-M Strehl','Location','NorthEastOutside');

figure;
hold on;
bar(1,mean(aicSpatFilterLMS),'FaceColor',[0 0 1]);
bar(2,mean(aicSpatFilterLminusM),'FaceColor',[1 0.5 0]);
bar(3,mean(aicStrehlLminusM),'FaceColor',[0.7 0.45 0]);
bar(4,mean(aicStrehlLMS));
bar(5,mean(aicDeltaPass2),'FaceColor',[0.7 0.7 0.7]);
bar(6,mean(aicSpatFilterLM),'FaceColor',[1 1 1]);
bar(7,mean(aicStrehlLM),'FaceColor',[1 1 1]);
set(gca,'XTick',[]);
set(gca,'FontSize',15);
ylabel('AIC');
legend('L+M-S','L-M','L-M Strehl','L+M-S Strehl','L+M 2cpd','L+M','Location','NorthEastOutside');

%% FIGURE 4C WITH EXTRA COMPARISON

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];

load([foldername 'wvMeanAndPredLM.mat']);
aicSpatFilterLM = aicAll;
load([foldername 'wvMeanAndPredLminusM.mat']);
aicSpatFilterLminusM = aicAll;
load([foldername 'wvMeanAndPredDonutx2.mat']);
aicSpatFilterLMS = aicAll;

figure; 
hold on;
boxplot([(aicSpatFilterLM-aicSpatFilterLMS)' ...
          (aicSpatFilterLM-aicSpatFilterLminusM)' ... 
          (aicSpatFilterLminusM-aicSpatFilterLMS)']);
ylim([-25 120]);
set(gca,'YTick',[-23.0259 0 23.0259 46.0517 69.0776 92.1034 115.1293]);
set(gca,'YTickLabel',{'10^-10' '1' '10^10' '10^20' '10^30' '10^40' '10^50'});
set(gca,'XTick',1:3);
set(gca,'XTickLabel',{'Blue-yellow' 'Red-green' 'Luminance'});

%% FIGURE 4C

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end
foldername = [dataPath 'data' slash 'PresavedFigureData' slash];

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
ylim([-25 120]);
set(gca,'YTick',[-23.0259 0 23.0259 46.0517 69.0776 92.1034 115.1293]);
set(gca,'YTickLabel',{'10^-10' '1' '10^10' '10^20' '10^30' '10^40' '10^50'});
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'Blue-yellow vs Luminance' 'Red-green vs Luminance'});
