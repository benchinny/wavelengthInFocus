% %%
% 
% subjNumAll = [1 3 5 10 16 17 18 20];
% 
% aicDeltaPass2 = [];
% 
% for i = 1:length(subjNumAll)
%     aic = ARCtestWvInFocusMeanZdeltaPassPlotStack(subjNumAll(i),'Lum',2);
%     aicDeltaPass2(i) = aic;
% end
% 
% save('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicDeltaPass2.mat','aicDeltaPass2');
% 
% %%
% 
% subjNumAll = [1 3 5 10 16 17 18 20];
% 
% aicDeltaPass4 = [];
% 
% for i = 1:length(subjNumAll)
%     aic = ARCtestWvInFocusMeanZdeltaPassPlotStack(subjNumAll(i),'Lum',4);
%     aicDeltaPass4(i) = aic;
% end
% 
% save('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicDeltaPass4.mat','aicDeltaPass4');

% %%
% 
% subjNumAll = [1 3 5 10 16 17 18 20];
% 
% aicStrehlLMS = [];
% 
% for i = 1:length(subjNumAll)
%     aic = ARCtestWvInFocusMeanZstrehlLMSeffectPlotStack(subjNumAll(i),'LMS');
%     aicStrehlLMS(i) = aic;
% end
% 
% save('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicStrehlLMS.mat','aicStrehlLMS');

% %%
% 
% subjNumAll = [1 3 5 10 16 17 18 20];
% 
% aicSpatFilterLMS = [];
% 
% for i = 1:length(subjNumAll)
%     aic = ARCtestWvInFocusMeanZspatFilterLMSeffectPlotStack(subjNumAll(i),'LMS');
%     aicSpatFilterLMS(i) = aic;
% end
% 
% save('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicSpatFilterLMS.mat','aicSpatFilterLMS');

% %%
% 
% subjNumAll = [1 3 5 10 16 17 18 20];
% 
% aicSpatFilterLM = [];
% 
% for i = 1:length(subjNumAll)
%     aic = ARCtestWvInFocusMeanZspatFilterLMSeffectPlotStack(subjNumAll(i),'LM');
%     aicSpatFilterLM(i) = aic;
% end
% 
% save('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicSpatFilterLM.mat','aicSpatFilterLM');
% 
%%

subjNumAll = [1 3 5 10 16 17 18 20];

aicStrehlLM = [];

for i = 1:length(subjNumAll)
    aic = ARCtestWvInFocusMeanZstrehlLMSeffectPlotStack(subjNumAll(i),'LM');
    aicStrehlLM(i) = aic;
end

save('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicStrehlLM.mat','aicStrehlLM');

% %%
% 
% subjNumAll = [1 3 5 10 16 17 18 20];
% 
% aicSpatFilterLminusM = [];
% 
% for i = 1:length(subjNumAll)
%     aic = ARCtestWvInFocusMeanZspatFilterLMSeffectPlotStack(subjNumAll(i),'LminusM');
%     aicSpatFilterLminusM(i) = aic;
% end
% 
% save('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicSpatFilterLminusM.mat','aicSpatFilterLminusM');

% %%
% 
% subjNumAll = [1 3 5 10 16 17 18 20];
% 
% aicStrehlLminusM = [];
% 
% for i = 1:length(subjNumAll)
%     aic = ARCtestWvInFocusMeanZstrehlLMSeffectPlotStack(subjNumAll(i),'LminusM');
%     aicStrehlLminusM(i) = aic;
% end
% 
% save('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicStrehlLminusM.mat','aicStrehlLminusM');

%%

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicDeltaPass2.mat')
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicDeltaPass4.mat')
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicSpatFilterLM.mat')
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicSpatFilterLminusM.mat')
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicSpatFilterLMS.mat')
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicStrehlLMS.mat')
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicStrehlLminusM.mat')
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/AICmodelComparisons/aicStrehlLM.mat')

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

%%

figure; 
hold on;
boxplot([(aicSpatFilterLM-aicSpatFilterLMS)' ...
          (aicSpatFilterLM-aicSpatFilterLminusM)' ... 
          (aicSpatFilterLminusM-aicSpatFilterLMS)']);
ylim([-25 95]);
set(gca,'YTick',[-23.0259 0 23.0259 46.0517 69.0776 92.1034]);
set(gca,'YTickLabel',{'10^-10' '1' '10^10' '10^20' '10^30' '10^40'});

%%

figure; 
hold on;
boxplot((aicSpatFilterLM-min([aicSpatFilterLMS; aicSpatFilterLminusM]))');
ylim([-25 95]);
set(gca,'YTick',[-23.0259 0 23.0259 46.0517 69.0776 92.1034]);
set(gca,'YTickLabel',{'10^-10' '1' '10^10' '10^20' '10^30' '10^40'});
