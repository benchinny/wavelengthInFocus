%%

subjNumAll = [1 3 5 10 16 17 18 20];

defocusAt550all = [];
rgbLumNormCndAll = [];
optDistCndAll = [];
subjNumVector = [];

for i = 1:length(subjNumAll)
    [wvInFocusCell, defocusAt550cell, optDistCnd, rgbLumNormCnd] = ARCnlz_lumEffectCheck(subjNumAll(i)+10);
    for j = 1:length(defocusAt550cell)
        defocusAt550all = [defocusAt550all; defocusAt550cell{j}];
        rgbLumNormCndAll = [rgbLumNormCndAll; repmat(rgbLumNormCnd(j,:),[length(defocusAt550cell{j}) 1])];
        optDistCndAll = [optDistCndAll; repmat(optDistCnd(j),[length(defocusAt550cell{j}) 1])];
        subjNumVector = [subjNumVector; repmat(subjNumAll(i),[length(defocusAt550cell{j}) 1])];
    end
end

%%

clear all; 

%%

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/lumControl/lumControlAnalysisData.mat');

%%

subjNumUnq = unique(subjNumVector);
rgbLumNormCndUnq = unique(rgbLumNormCndAll,'rows');
defocusAt550mean = [];

for i = 1:length(subjNumUnq)
    for j = 1:size(rgbLumNormCndUnq)
        ind = abs(rgbLumNormCndAll(:,1)-rgbLumNormCndUnq(j,1))<0.001 & ...
              abs(rgbLumNormCndAll(:,2)-rgbLumNormCndUnq(j,2))<0.001 & ...
              abs(rgbLumNormCndAll(:,3)-rgbLumNormCndUnq(j,3))<0.001 & ...
              subjNumVector==subjNumUnq(i);
        defocusAt550mean(i,j) = mean(defocusAt550all(ind));
    end
end

defocusAt550diff = defocusAt550mean(:,2)-defocusAt550mean(:,1);
defocusAt550meanAcrossSubj = mean(defocusAt550mean);
defocusAt550stdAcrossSubj = std(defocusAt550mean)/sqrt(length(subjNumUnq));

figure;
hold on;
bar([0.5 0.8],-defocusAt550meanAcrossSubj,'w','LineWidth',1);
errorbar([0.5 0.8],-defocusAt550meanAcrossSubj,defocusAt550stdAcrossSubj,'k.','LineWidth',1);
axis square;
xlim([0.3 1]);
ylim([0 3]);
set(gca,'FontSize',15);
set(gca,'XTick',[0.5 0.8 1.0]);
set(gca,'XTickLabel',{'0.5' '0.8' '1.0'});
xlabel('Luminance (cd/m^2)');
ylabel('Defocus at 550nm (D)');

%%

rgbLumNormCndUnq = unique(rgbLumNormCndAll,'rows');
defocusAt550mean = [];

for i = 1:size(rgbLumNormCndUnq,1)
    ind = abs(rgbLumNormCndAll(:,1)-rgbLumNormCndUnq(i,1))<0.001 & ...
          abs(rgbLumNormCndAll(:,2)-rgbLumNormCndUnq(i,2))<0.001 & ...
          abs(rgbLumNormCndAll(:,3)-rgbLumNormCndUnq(i,3))<0.001;
    defocusAt550mean(i) = mean(defocusAt550all(ind));
    defocusAt550std(i) = std(defocusAt550all(ind))/sqrt(sum(ind));
end

figure;
hold on;
bar([0.5 0.8],-defocusAt550mean,'w','LineWidth',1);
errorbar([0.5 0.8],-defocusAt550mean,defocusAt550std,'k.','LineWidth',1);
axis square;
xlim([0.3 1]);
ylim([0 3]);
set(gca,'FontSize',15);
set(gca,'XTick',[0.5 0.8 1.0]);
set(gca,'XTickLabel',{'0.5' '0.8' '1.0'});
xlabel('Luminance (cd/m^2)');
ylabel('Defocus at 550nm (D)');
