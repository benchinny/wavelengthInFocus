%% FIGURE 5C IN MANUSCRIPT

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

% ALL SUBJECT NUMBERS
subjNumAll = [1 3 5 10 16 17 18 20];

dprimeRatioAll = []; % D-PRIME RATIOS
for i = 1:length(subjNumAll)
    [~,~,~,dprime,~,~,~,~,~,~,~] = ARCacuityAnalyzeDataOnly(subjNumAll(i),0,dataPath);
    % RATIO OF BEST D-PRIME TO D-PRIME AT 0
    dprimeRatioAll(i) = max(dprime)/dprime(5);
end

figure;
boxplot(dprimeRatioAll);
set(gca,'FontSize',15);
hold on;
for i = 1:length(dprimeRatioAll)
    if abs(dprimeRatioAll(i)-1)<0.01
        plot(1+i*0.1,dprimeRatioAll(i),'ko','MarkerSize',10,'MarkerFaceColor',[0.56 0 1]);
    else
        plot(1,dprimeRatioAll(i),'ko','MarkerSize',10,'MarkerFaceColor',[0.56 0 1]);
    end
end
ylabel('d''_{max}/d''_{0}');
ylim([0.9 2.6]);
set(gca,'YTick',[1 1.5 2 2.5]);

% CALCULATE, FOR EACH SUBJECT, PROBABILITY OF THE PEAK BEING AT 0 (STATS
% REPORTED IN MANUSCRIPT)
pPeakAt0 = ARCnlz_fig5Cstats(dataPath);

% DISPLAY PROBABILITIES
for i = 1:length(pPeakAt0)
    display(['Subject S' num2str(subjNumAll(i)) ' probability of peaking at 0 = ' num2str(pPeakAt0(i),3)]);
end