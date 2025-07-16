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

% MAKE EXAMPLE CONE IMAGES FOR SUBJECT 2 AT INFORMATIVE WAVELENGTHS

wave = 380:4:780;
wv2vis = [620 540 460];
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

redblue = make_yellow_blue_colormap(1);
axisLims = [min([coneImgFilteredEgLM(:); coneImgFilteredEgLMS(:)]) ...
            max([coneImgFilteredEgLM(:); coneImgFilteredEgLMS(:)])];
stimCell = {'More blue' 'Purple' 'More red'};
bColorBar = false;

figure;
set(gcf,'Position',[29 90 890 908]);
for i = 1:length(wv2vis)
    subplot(3,3,(i-1)*3+3);
    imagesc(squeeze(coneImgFilteredEgLMS(:,:,1,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'LineWidth',2);
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    if i==1 && bColorBar
       colorbar;
    end
    clim(axisLims(2).*[-1 1])
    if bColorBar
        title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{1} ', LMS']);
    end
    subplot(3,3,(i-1)*3+2);
    imagesc(squeeze(absorptionsS_LMS(:,:,1,i)));
    axis square; 
    colormap(redblue); 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]); 
    set(gca,'LineWidth',2);
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);
    if bColorBar
       title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{1} ', S image']);
    end
    subplot(3,3,(i-1)*3+1);
    imagesc(squeeze(absorptionsLM_LMS(:,:,1,i)));
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]); 
    set(gca,'LineWidth',2);
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);
    if bColorBar
        title('L+M image');
    end
end

figure;
set(gcf,'Position',[29 90 890 908]);
for i = 1:length(wv2vis)
    subplot(3,3,(i-1)*3+3);
    imagesc(squeeze(coneImgFilteredEgLMS(:,:,3,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'LineWidth',2);
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    if i==1 && bColorBar
       colorbar;
    end
    clim(axisLims(2).*[-1 1])
    if bColorBar
        title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{3} ', LMS']);
    end
    subplot(3,3,(i-1)*3+2);
    imagesc(squeeze(absorptionsS_LMS(:,:,3,i)));
    axis square; 
    colormap(redblue); 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]); 
    set(gca,'LineWidth',2);
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);  
    if bColorBar
        title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{3} ', S image']);
    end
    subplot(3,3,(i-1)*3+1);
    imagesc(squeeze(absorptionsLM_LMS(:,:,3,i)));
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]); 
    set(gca,'LineWidth',2);
    colormap(redblue); 
    xlim([50 120]); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);
    if bColorBar
       title('L+M image');
    end
end

%% FIRST FINCH PLOT

dataFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/';
FinchData = readtable([dataFolder 'Finch_data_RedVsGreenErrorBars' '.csv']);
lumRatio = 0:0.125:1;
accommodation = [0.179 0.144 0.106 0.055 0.029 0.025 -0.011 -0.072 -0.175];
ubAcc = [0.204 0.172 0.135 0.0844 0.057 0.0532 0.0128 -0.0443 -0.149];
lbAcc = [0.149 0.120 0.075 0.0285 -0.0028 -0.0046 -0.0367 -0.0979 -0.202];
figure;
errorbar(lumRatio,accommodation,ubAcc-accommodation,accommodation-lbAcc,'ko','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
hold on;
plot(xlim,0.273.*[1 1],'r-','LineWidth',1);
plot(xlim,-0.335.*[1 1],'g-','LineWidth',1);
axis square;
set(gca,'FontSize',15);
xlabel('Green/Red Luminance');
ylabel('Accommodation');
ylim([-1 0.3]);

%% SECOND FINCH PLOT

dataFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/';
FinchData = readtable([dataFolder 'Finch_data_RedVsBlueErrorBars' '.csv']);
lumRatio = 0:0.125:1;
accommodation = [0.178 0.223 0.202 0.150 0.0890 0.031 -0.061 -0.209 -0.352];
ubAcc = [0.204 0.251 0.231 0.178 0.124 0.064 -0.034 -0.177 -0.316];
lbAcc = [0.146 0.191 0.163 0.122 0.052 -0.005 -0.089 -0.236 -0.387];
figure;
errorbar(lumRatio,accommodation,ubAcc-accommodation,accommodation-lbAcc,'ko','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
hold on;
plot(xlim,0.276.*[1 1],'r-','LineWidth',1);
plot(xlim,-0.873.*[1 1],'b-','LineWidth',1);
axis square;
set(gca,'FontSize',15);
xlabel('Blue/Red Luminance');
ylabel('Accommodation');
ylim([-1 0.3]);

