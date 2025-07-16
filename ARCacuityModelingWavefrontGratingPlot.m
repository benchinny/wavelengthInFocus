%%

clear;

%%

folderPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_Dec_18/';

subjNum = 9;

if subjNum==10
filenames = { ...
             'dprimeRedS10.mat' ...
             'dprimeGreenS10.mat' ...
             'dprimeBlueS10.mat' ...           
             };
elseif subjNum==3
filenames = { ...
             'dprimeRedS3.mat' ...
             'dprimeGreenS3.mat' ...
             'dprimeBlueS3.mat' ...           
             }; 
elseif subjNum==9
filenames = { ...
             'dprimeRedS9.mat' ...
             'dprimeGreenS9.mat' ...
             'dprimeBlueS9.mat' ...           
             }; 
end

[defocusLCAmeasured, q1best, q2best, q3best,defocusLCAmeasuredBoots,Dgreen,dprimeFitAll,PCfitSupport] = ARCacuAnalysisLCA(subjNum,0,0);

figure;
hold on;
colorLetters = 'rgb';
for i = 1:3
    load([folderPath filenames{i}]);
    scaleFactor = max(dprimeFitAll(i,:))/max(dprimeMetric);
    plot(PCfitSupport+2.5,dprimeFitAll(i,:)./scaleFactor,[colorLetters(i) '-']);
    plot(fliplr(defocusForStim+defocusOrigScaled),dprimeMetric,[colorLetters(i) '--']);
end
axis square;
set(gca,'FontSize',12);
xlabel('Stimulus distance');
ylabel('D-prime');

