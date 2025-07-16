%% COMPARING MEAN TRACES

sn = [18 23 24 26 27 28 29 30 32 33];
frmDuration = 0.033;
% lumCutoff = [0.2 0.65; 0.65 1.1; 1.1 1.75];
lumCutoff = [0 5];

for j = 1:size(lumCutoff,1)
    figure;
    set(gcf,'Position',[105 245 1479 647]);    
    for i = 1:length(sn) 
        subplot(2,5,i);
        set(gca,'FontSize',15);
        hold on;
        [x3stackRed,x3stackBlue,x3stackMixed,x3stackMixedMoreRed,x3stackMixedMoreBlue,meanRedPerTrial,meanBluePerTrial,meanMixedPerTrial,meanMixedMoreRedPerTrial,meanMixedMoreBluePerTrial,tInterp,AFCp,indCndCell] = ARCnlz_redVsBlueVsMixed(sn(i),0,lumCutoff(j,:));
        meanAnchor = mean(x3stackRed(:,1));
        tSamples = [0:frmDuration:(size(x3stackMixed,2)-1)*frmDuration];
        xStackCIred = quantile(x3stackRed,[0.16 0.84])-meanAnchor;
        xStackCIblue = quantile(x3stackBlue,[0.16 0.84])-meanAnchor;
        xStackCImixed = quantile(x3stackMixed,[0.16 0.84])-meanAnchor;
        fill([tSamples fliplr(tSamples)],[xStackCIred(1,:) fliplr(xStackCIred(2,:))],[1 0 0],'FaceAlpha',0.1,'EdgeColor','none');
        fill([tSamples fliplr(tSamples)],[xStackCIblue(1,:) fliplr(xStackCIblue(2,:))],[0 0 1],'FaceAlpha',0.1,'EdgeColor','none');
        fill([tSamples fliplr(tSamples)],[xStackCImixed(1,:) fliplr(xStackCImixed(2,:))],[1 0 1],'FaceAlpha',0.1,'EdgeColor','none');
        plot([0:frmDuration:(size(x3stackRed,2)-1)*frmDuration],mean(x3stackRed)-meanAnchor,'-','Color',[1 0 0],'LineWidth',2);
        plot([0:frmDuration:(size(x3stackBlue,2)-1)*frmDuration],mean(x3stackBlue)-meanAnchor,'-','Color',[0 0 1],'LineWidth',2);
        plot([0:frmDuration:(size(x3stackMixed,2)-1)*frmDuration],mean(x3stackMixed)-meanAnchor,'-','Color',[1 0 1],'LineWidth',2);
        display(['Done with subject: ' num2str(i)]);
        axis square;
        xlim([0 3]);
        ylim([-3 3]);
        xlabel('Time (s)'); ylabel('Relative Power (Diopters)'); title(['Subject ' num2str(sn(i))]);
    end
end

%% HISTOGRAMS OF MEANS FOR EACH TRIAL

sn = [18 23 24 26 27 28 29 30 32 33];
frmDuration = 0.033;
lumCutoff = [0.2 1.75];
meanDiffPureRedVsBlue = [];

for j = 1:size(lumCutoff,1)
    figure;
    set(gcf,'Position',[105 245 1479 647]);    
    for i = 1:length(sn) 
        subplot(2,5,i);
        set(gca,'FontSize',15);
        hold on;
        [x3stackRed,x3stackBlue,x3stackMixed,x3stackMixedMoreRed,x3stackMixedMoreBlue,meanRedPerTrial,meanBluePerTrial,meanMixedPerTrial,meanMixedMoreRedPerTrial,meanMixedMoreBluePerTrial,tInterp,AFCp,indCndCell] = ARCnlz_redVsBlueVsMixed(sn(i),0,lumCutoff(j,:));
        histogram(meanMixedPerTrial-mean(meanRedPerTrial),'FaceColor',[1 0 1]);
        plot((mean(meanRedPerTrial)-mean(meanRedPerTrial)).*[1 1],ylim,'r-','LineWidth',1.5);
        plot((mean(meanBluePerTrial)-mean(meanRedPerTrial)).*[1 1],ylim,'b-','LineWidth',1.5);
        display(['Done with subject: ' num2str(i)]);
        axis square;
        xlabel('Measurement'); ylabel('Count'); 
        meanDiffPureRedVsBlue(i) = mean(meanRedPerTrial)-mean(meanBluePerTrial);
    end
end

%% COMPARING MEAN TRACES FOR MORE RED VS MORE BLUE

sn = [18 23 24 26 27 29 32 33];
frmDuration = 0.033;

figure;
set(gcf,'Position',[105 245 1479 647]);    
for i = 1:length(sn) 
    subplot(2,4,i);
    set(gca,'FontSize',15);
    hold on;
    [x3stackRed,x3stackBlue,x3stackMixed,x3stackMixedMoreRed,x3stackMixedMoreBlue,meanRedPerTrial,meanBluePerTrial,meanMixedPerTrial,meanMixedMoreRedPerTrial,meanMixedMoreBluePerTrial,tInterp,AFCp,indCndCell] = ARCnlz_redVsBlueVsMixed(sn(i),0,[0 2]);
    meanAnchor = mean(x3stackRed(:,1));
    tSamples = [0:frmDuration:(size(x3stackMixed,2)-1)*frmDuration];
    xStackCIred = quantile(x3stackRed,[0.16 0.84])-meanAnchor;
    xStackCIblue = quantile(x3stackBlue,[0.16 0.84])-meanAnchor;
    xStackCImixedMoreRed = quantile(x3stackMixedMoreRed,[0.16 0.84])-meanAnchor;
    xStackCImixedMoreBlue = quantile(x3stackMixedMoreBlue,[0.16 0.84])-meanAnchor;
    fill([tSamples fliplr(tSamples)],[xStackCIred(1,:) fliplr(xStackCIred(2,:))],[1 0 0],'FaceAlpha',0.1,'EdgeColor','none');
    fill([tSamples fliplr(tSamples)],[xStackCIblue(1,:) fliplr(xStackCIblue(2,:))],[0 0 1],'FaceAlpha',0.1,'EdgeColor','none');
    fill([tSamples fliplr(tSamples)],[xStackCImixedMoreRed(1,:) fliplr(xStackCImixedMoreRed(2,:))],[1 0 0.5],'FaceAlpha',0.1,'EdgeColor','none');
    fill([tSamples fliplr(tSamples)],[xStackCImixedMoreBlue(1,:) fliplr(xStackCImixedMoreBlue(2,:))],[0.5 0 1],'FaceAlpha',0.1,'EdgeColor','none');
    plot([0:frmDuration:(size(x3stackRed,2)-1)*frmDuration],mean(x3stackRed)-meanAnchor,'-','Color',[1 0 0],'LineWidth',2);
    plot([0:frmDuration:(size(x3stackBlue,2)-1)*frmDuration],mean(x3stackBlue)-meanAnchor,'-','Color',[0 0 1],'LineWidth',2);
    plot([0:frmDuration:(size(x3stackMixedMoreRed,2)-1)*frmDuration],mean(x3stackMixedMoreRed)-meanAnchor,'-','Color',[1 0 0.5],'LineWidth',2);
    plot([0:frmDuration:(size(x3stackMixedMoreBlue,2)-1)*frmDuration],mean(x3stackMixedMoreBlue)-meanAnchor,'-','Color',[0.5 0 1],'LineWidth',2);
    display(['Done with subject: ' num2str(i)]);
    axis square;
    xlim([0 3]);
    ylim([-3 3]);
    xlabel('Time (s)'); ylabel('Relative Power (Diopters)');     
end

%% HISTOGRAMS OF MEANS FOR EACH TRIAL FOR MIXED MORE RED VS MORE BLUE

sn = [18 23 24 26 27 29 32 33];
frmDuration = 0.033;
lumCutoff = [0.2 1.75];

for j = 1:size(lumCutoff,1)
    figure;
    set(gcf,'Position',[105 245 1479 647]);    
    for i = 1:length(sn) 
        subplot(2,4,i);
        set(gca,'FontSize',15);
        hold on;
        [x3stackRed,x3stackBlue,x3stackMixed,x3stackMixedMoreRed,x3stackMixedMoreBlue,meanRedPerTrial,meanBluePerTrial,meanMixedPerTrial,meanMixedMoreRedPerTrial,meanMixedMoreBluePerTrial,tInterp,AFCp,indCndCell] = ARCnlz_redVsBlueVsMixed(sn(i),0,lumCutoff(j,:));
        histogram(meanMixedMoreRedPerTrial-mean(meanRedPerTrial),'FaceColor',[1 0 0.5]);
        histogram(meanMixedMoreBluePerTrial-mean(meanRedPerTrial),'FaceColor',[0.5 0 1]);
        plot((mean(meanRedPerTrial)-mean(meanRedPerTrial)).*[1 1],ylim,'r-','LineWidth',1.5);
        plot((mean(meanBluePerTrial)-mean(meanRedPerTrial)).*[1 1],ylim,'b-','LineWidth',1.5);
        display(['Done with subject: ' num2str(i)]);
        axis square;
        xlabel('Measurement'); ylabel('Count');     
    end
end

%% COMPARING MEAN TRACES

sn = [11 16];
frmDuration = 0.033;
lumCutoff = [0 5];

for j = 1:size(lumCutoff,1)
    figure;
    set(gcf,'Position',[105 245 1479 647]);    
    for i = 1:length(sn) 
        subplot(2,4,i);
        set(gca,'FontSize',15);
        hold on;
        [x3stackRed,x3stackBlue,x3stackMixed,x3stackMixedMoreRed,x3stackMixedMoreBlue,meanRedPerTrial,meanBluePerTrial,meanMixedPerTrial,meanMixedMoreRedPerTrial,meanMixedMoreBluePerTrial,tInterp,AFCp,indCndCell] = ARCnlz_redVsBlueVsMixed(sn(i),0,lumCutoff(j,:));
        meanAnchor = mean(x3stackRed(:,1));
        tSamples = [0:frmDuration:(size(x3stackMixed,2)-1)*frmDuration];
        xStackCIred = quantile(x3stackRed,[0.16 0.84])-meanAnchor;
        xStackCIblue = quantile(x3stackBlue,[0.16 0.84])-meanAnchor;
        xStackCImixed = quantile(x3stackMixed,[0.16 0.84])-meanAnchor;
        fill([tSamples fliplr(tSamples)],[xStackCIred(1,:) fliplr(xStackCIred(2,:))],[1 0 0],'FaceAlpha',0.1,'EdgeColor','none');
        fill([tSamples fliplr(tSamples)],[xStackCIblue(1,:) fliplr(xStackCIblue(2,:))],[0 0 1],'FaceAlpha',0.1,'EdgeColor','none');
        fill([tSamples fliplr(tSamples)],[xStackCImixed(1,:) fliplr(xStackCImixed(2,:))],[1 0 1],'FaceAlpha',0.1,'EdgeColor','none');
        plot([0:frmDuration:(size(x3stackRed,2)-1)*frmDuration],mean(x3stackRed)-meanAnchor,'-','Color',[1 0 0],'LineWidth',2);
        plot([0:frmDuration:(size(x3stackBlue,2)-1)*frmDuration],mean(x3stackBlue)-meanAnchor,'-','Color',[0 0 1],'LineWidth',2);
        plot([0:frmDuration:(size(x3stackMixed,2)-1)*frmDuration],mean(x3stackMixed)-meanAnchor,'-','Color',[1 0 1],'LineWidth',2);
        display(['Done with subject: ' num2str(i)]);
        axis square;
        xlim([0 3]);
        ylim([-3 3]);
        xlabel('Time (s)'); ylabel('Relative Power (Diopters)');     
    end
end


%%

clear;