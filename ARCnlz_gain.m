%% COMPARE WEIGHTING MODEL AND SWITCHING MODEL

sn = [18 23 24 26 27 28 29 30 32 33];
lumBins = [0 0.3; 0.4 0.7];

weightsSall = [];
trialCount = [];

for i = 1:length(sn)
    for j = 1:size(lumBins,1)
        [weightsS1_x, weightsRBS1_y, rhoFull, rhoNoColor, rhoColor, aicLin, aicNoColor, weightsRBSci, trialMeans, meanChange, deltaR, deltaB, deltaS] = ARCnlz_linearModelGainOnly(sn(i),false,1000,lumBins(j,:));
        weightsSall(i,j) = weightsS1_x;
        trialCount(i,j) = length(meanChange)
        display(['Done with subject ' num2str(i)]);
    end
end

%%

figure; 
plot(weightsSall','o','LineWidth',1.5,'MarkerSize',15,'MarkerFaceColor','w');
axis square;
ylim([0 1.5]);
xlim([0.5 2.5]);
set(gca,'FontSize',15);
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'Lower luminance' 'Higher luminance'});
xlabel('Condition');
ylabel('Gain');
