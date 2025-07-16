%% RECORD OF BEST FIT CONE WEIGHTS

LMSweights = [0.5467    0.2533   -1.0000; ...
              0.4725    0.5775   -1.0000; ...
              0.7650    0.1350   -0.5000; ...
              1.1475    0.2025   -1.0000; ...
              0.2267    0.5733   -1.0000; ...
              0.3008    0.6492   -1.0000; ...
              1.1050    0.1950   -1.0000; ...
              1.1433    0.2567   -0.2500];

LMweights = [0.8250    0.1750         0; ...
             0.8500    0.1500         0; ...
             0.8500    0.1500         0; ...
             0.8500    0.1500         0; ...
             0.2500    0.7500         0; ...
             0.8250    0.1750         0; ...
             0.8500    0.1500         0; ...
             0.7500    0.2500         0];

LminusMweights = [0.3250   -0.1750         0; ...
                  0.3750   -0.1250         0; ...
                  0.4125   -0.0875         0; ...
                  0.8000   -0.2000         0; ...
                  0.7000   -0.3000         0; ...
                  0.4000   -0.1000         0; ...
                  0.6250   -0.3750         0; ...
                  0.4875   -0.0125         0];

%% GET WAVELENGTH-IN-FOCUS PREDICTIONS FROM EACH MODEL

subjNumAll = [1 3 5 10 16 17 18 20];
stimNum = [1 3 8 7 6];
wvInFocusLMS = [];
wvInFocusLM = [];
wvInFocusLminusM = [];
peakCorrLM = [];
peakCorrLMS = [];
peakCorrLminusM = [];

for i = 1:size(LMSweights,1)
    for j = 1:length(stimNum)
        [wvInFocus, coneImgFilteredEg, coneImgOrigFilteredEg, wvInFocus2, wave, peakCorr] = ...
         ARCwvInFocusConesMeanZsandbox(subjNumAll(i),stimNum(j),LMSweights(i,:),[]);
        wvInFocusLMS(i,j,:) = wvInFocus;
        peakCorrLMS(i,j,:) = peakCorr;
        [wvInFocus, coneImgFilteredEg, coneImgOrigFilteredEg, wvInFocus2, wave, peakCorr] = ...
         ARCwvInFocusConesMeanZsandbox(subjNumAll(i),stimNum(j),LMweights(i,:),[]);
        wvInFocusLM(i,j,:) = [wvInFocus wvInFocus2];
        peakCorrLM(i,j,:) = peakCorr;
        [wvInFocus, coneImgFilteredEg, coneImgOrigFilteredEg, wvInFocus2, wave, peakCorr] = ...
         ARCwvInFocusConesMeanZsandbox(subjNumAll(i),stimNum(j),LminusMweights(i,:),[]);
        wvInFocusLminusM(i,j,:) = wvInFocus;
        peakCorrLminusM(i,j,:) = peakCorr;
        display(['Weight index ' num2str(i) 'Stimulus index ' num2str(j)]);
    end
    disp(['Subj ' num2str(i) ' done']);
end

%% SAVE RESULTS OF PREVIOUS BLOCK

save('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/helperFiles/wvInFocusAllStimAndSubj.mat', ...
     'wvInFocusLMS','wvInFocusLminusM','wvInFocusLM','peakCorrLMS','peakCorrLminusM','peakCorrLM');

%% MAKE EXAMPLE CONE IMAGES FOR SUBJECT 2 AT INFORMATIVE WAVELENGTHS

wave = 380:4:780;
wv2vis = [620 540 460 420];
stimInd = [1 8 6];
coneImgFilteredEgLMS = [];
coneImgFilteredEgLM = [];
absorptionsS_LMS = [];
absorptionsLM_LMS = [];

for i = 1:length(stimInd)
    for j = 1:length(wv2vis)
        [wvInFocus, coneImgFilteredEg, coneImgOrigFilteredEg, wvInFocus2, wave, peakCorr, absorptions] = ...
         ARCwvInFocusConesMeanZsandbox(3,stimInd(i),LMSweights(2,:),find(wave==wv2vis(j)));
         coneImgFilteredEgLMS(:,:,i,j) = coneImgFilteredEg;
         absorptionsS_LMS(:,:,i,j) = absorptions(:,:,3);
         absorptionsLM_LMS(:,:,i,j) = absorptions(:,:,1)+absorptions(:,:,2);

        [wvInFocus, coneImgFilteredEg, coneImgOrigFilteredEg, wvInFocus2, wave, peakCorr, absorptions] = ...
         ARCwvInFocusConesMeanZsandbox(3,stimInd(i),[LMSweights(2,1) LMSweights(2,2) 0],find(wave==wv2vis(j))); 
        coneImgFilteredEgLM(:,:,i,j) = coneImgFilteredEg;
    end
end

%% PLOT DIFFERENT CONE IMAGES FOR BUILDING MODEL INTUITIONS

redblue = make_red_blue_colormap(1);
axisLims = [min([coneImgFilteredEgLM(:); coneImgFilteredEgLMS(:)]) ...
            max([coneImgFilteredEgLM(:); coneImgFilteredEgLMS(:)])];
stimCell = {'More blue' 'Purple' 'More red'};

figure;
set(gcf,'Position',[16 103 558 868]);
for i = 1:length(wv2vis)
    subplot(4,2,(i-1)*2+1);
    imagesc(squeeze(coneImgFilteredEgLMS(:,:,1,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    colorbar; 
    clim(axisLims(2).*[-1 1])
    title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{1} ', LMS']);
    subplot(4,2,(i-1)*2+2);
    imagesc(squeeze(coneImgFilteredEgLM(:,:,1,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1])    
    title('LM');
end

figure;
set(gcf,'Position',[551 99 558 868]);
for i = 1:length(wv2vis)
    subplot(4,2,(i-1)*2+1);
    imagesc(squeeze(coneImgFilteredEgLMS(:,:,2,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1])
    title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{2} ', LMS']);
    subplot(4,2,(i-1)*2+2);
    imagesc(squeeze(coneImgFilteredEgLM(:,:,2,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1])   
    title('LM');
end

figure;
set(gcf,'Position',[1073 96 558 868]);
for i = 1:length(wv2vis)
    subplot(4,2,(i-1)*2+1);
    imagesc(squeeze(coneImgFilteredEgLMS(:,:,3,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1])
    title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{3}]);
    subplot(4,2,(i-1)*2+2);
    imagesc(squeeze(coneImgFilteredEgLM(:,:,3,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);
end

figure;
set(gcf,'Position',[16 103 558 868]);
for i = 1:length(wv2vis)
    subplot(4,2,(i-1)*2+1);
    imagesc(squeeze(absorptionsS_LMS(:,:,1,i)));
    axis square; 
    colormap(redblue); 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);  
    title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{1} ', S image']);
    subplot(4,2,(i-1)*2+2);
    imagesc(squeeze(absorptionsLM_LMS(:,:,1,i)));
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);
    title('L+M image');
end

figure;
set(gcf,'Position',[551 99 558 868]);
for i = 1:length(wv2vis)
    subplot(4,2,(i-1)*2+1);
    imagesc(squeeze(absorptionsS_LMS(:,:,2,i)));
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);  
    title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{2} ', S image']);
    subplot(4,2,(i-1)*2+2);
    imagesc(squeeze(absorptionsLM_LMS(:,:,2,i)));
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);
    title('L+M image');
end

figure;
set(gcf,'Position',[1073 96 558 868]);
for i = 1:length(wv2vis)
    subplot(4,2,(i-1)*2+1);
    imagesc(squeeze(absorptionsS_LMS(:,:,3,i)));
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);  
    title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{3} ', S image']);
    subplot(4,2,(i-1)*2+2);
    imagesc(squeeze(absorptionsLM_LMS(:,:,3,i)));
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);    
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);
    title('L+M image');
end
%% PLOT WAVELENGTH-IN-FOCUS FOR BLUE-YELLOW OPPONENT MODEL

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/helperFiles/wvInFocusAllModel.mat','wvInFocusAllTable','peakCorr8LMS','peakCorr8LM','peakCorr1LMS','peakCorr1LM','peakCorr6LMS','peakCorr6LM');

figure;
hold on;
for i = 1:8
    wvInFocusBlue = wvInFocusAllTable.("L+M-S1");
    wvInFocusPurple = wvInFocusAllTable.("L+M-S8");
    wvInFocusRed = wvInFocusAllTable.("L+M-S6");
    plot([wvInFocusBlue(i) wvInFocusPurple(i) wvInFocusRed(i)],'LineWidth',1);
end
plot(mean([wvInFocusBlue wvInFocusPurple wvInFocusRed]),'k-','LineWidth',2);
plot([0.5 3.5],[620 620],'r--','LineWidth',1);
plot([0.5 3.5],[460 460],'b--','LineWidth',1);
set(gca,'FontSize',15);
ylim([430 650]);
xlim([0.75 3.25]);
set(gca,'XTick',[1 2 3]);
set(gca,'XTickLabel',{'Bluest' 'Purple' 'Reddest'});
ylabel('Wavelength in focus (nm)');
legend({'S1' 'S2' 'S3' 'S4' 'S5' 'S6' 'S7' 'S8' 'Average'},'Location','SouthEast');
title('No green model predictions');

%% PLOT PEAK CORRELATION AS A FUNCTION OF WAVELENGTH FOR ALL SUBJECTS

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/helperFiles/wvInFocusAllModel.mat','wvInFocusAllTable','peakCorr8LMS','peakCorr8LM','peakCorr1LMS','peakCorr1LM','peakCorr6LMS','peakCorr6LM');
wave = 380:4:780;

figure;
set(gcf,'Position',[263 365 1287 574]);
for i = 1:8
   subplot(2,4,i);
   plot(wave,peakCorr1LMS(i,:),'b-','LineWidth',1); 
   hold on; 
   plot(wave,peakCorr1LM(i,:),'k-','LineWidth',1)
   if i==1
       xlabel('Wavelength in focus (nm)');
       ylabel('Normalized cross-correlation');
   end   
   set(gca,'FontSize',15);
   axis square;
   title('Bluest');
end

figure;
set(gcf,'Position',[263 365 1287 574]);
for i = 1:8
   subplot(2,4,i);
   plot(wave,peakCorr6LMS(i,:),'b-','LineWidth',1); 
   hold on; 
   plot(wave,peakCorr6LM(i,:),'k-','LineWidth',1)
   if i==1
       xlabel('Wavelength in focus (nm)');
       ylabel('Normalized cross-correlation');
   end
   set(gca,'FontSize',15);
   axis square;   
   title('Reddest');
end

figure;
set(gcf,'Position',[263 365 1287 574]);
for i = 1:8
   subplot(2,4,i);
   plot(wave,peakCorr8LMS(i,:),'b-','LineWidth',1); 
   hold on; 
   plot(wave,peakCorr8LM(i,:),'k-','LineWidth',1)
   if i==1
       xlabel('Wavelength in focus (nm)');
       ylabel('Normalized cross-correlation');
   end
   set(gca,'FontSize',15);
   axis square; 
   title('Purple');
end

%% PLOT CORRELATION AS A FUNCTION OF WAVELENGTH FOR DIFFERENT STIMULI

wave = 380:4:780;
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/helperFiles/wvInFocusAllStimAndSubj.mat', ...
'wvInFocusLMS','wvInFocusLminusM','wvInFocusLM','peakCorrLMS','peakCorrLminusM','peakCorrLM');

subjNumInd = 2;

lumRatio = [0.25 0 1.00; ...
            0.50 0 1.00; ...
            1.00 0 1.00; ...
            1.00 0 0.50; ...
            1.00 0 0.25];

figure; 
hold on;
for i = 1:size(peakCorrLM,2)
    plot(wave,squeeze(peakCorrLM(subjNumInd,i,:)),'LineWidth',1,'Color',lumRatio(i,:));
end
yLimTmp = ylim;
for i = 1:size(peakCorrLM,2)
    [peakVal,indPeak] = max(peakCorrLM(subjNumInd,i,:));
    wvInFocus = wave(indPeak);
    plot([wvInFocus wvInFocus],[yLimTmp(1) peakVal],'--','LineWidth',1,'Color',lumRatio(i,:));
end
xlim([400 700]);
set(gca,'FontSize',15);

figure; 
hold on;
for i = 1:size(peakCorrLMS,2)
    plot(wave,squeeze(peakCorrLMS(subjNumInd,i,:)),'LineWidth',1,'Color',lumRatio(i,:));
end
yLimTmp = ylim;
for i = 1:size(peakCorrLMS,2)
    [peakVal,indPeak] = max(peakCorrLMS(subjNumInd,i,:));
    wvInFocus = wave(indPeak);
    plot([wvInFocus wvInFocus],[yLimTmp(1) peakVal],'--','LineWidth',1,'Color',lumRatio(i,:));
end
xlim([400 700]);
set(gca,'FontSize',15);

% figure; 
% plot(wave,squeeze(peakCorrLminusM(1,:,:)));

