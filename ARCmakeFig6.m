%% SCRIPT FOR MAKING FIGURE 6 (ILLUSTRATE MODEL INTUITION)

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
subjNumInd = 2; % SUBJECT TO LOOK AT
% FOLDER NAME WITH PRESAVED DATA
foldername = fullfile(dataPath,'data','PresavedFigureData');

% LIST OF ALL SUBJECTS
subjNumAll = [1 3 5 10 16 17 18 20]; 

% LOAD WEIGHTS FROM PRESAVED DATA

% FOR L+M
load(fullfile(foldername,'wvMeanAndPredLM.mat'),'wLpropMinAll','wLMminAll');
wL = wLMminAll(2)*wLpropMinAll(2);
wM = wLMminAll(2)-wL;
LMweights = [wL wM 0];

% FOR BLUE-YELLOW
load(fullfile(foldername,'wvMeanAndPredLMS.mat'),'wLpropMinAll','wLMminAll','wSall');
wL = wLMminAll(2)*wLpropMinAll(2);
wM = wLMminAll(2)-wL;
LMSweights = [wL wM wSall(subjNumInd)];

% FOR RED-GREEN
load(fullfile(foldername,'wvMeanAndPredLminusM.mat'),'wLpropMinAll','wLMminAll');
wL = wLMminAll(2)*wLpropMinAll(2);
wM = wLMminAll(2)-wL;
LminusMweights = [wL wM 0];

%% MAKE EXAMPLE CONE IMAGES FOR SUBJECT 2 AT INFORMATIVE WAVELENGTHS

wave = 380:4:780; % WAVELENGTHS TO SAMPLE AT
wv2vis = [616 552 468]; % WAVELENGTHS TO VISUALIZE
stimInd = [1 8 6]; % STIMULUS INDICES CORRESPONDING TO STIMULI OF INTEREST
coneImgEgLMS = []; % FOR SAVING OUT CONE IMAGES FOR LMS
coneImgEgLM = []; % FOR SAVING OUT CONE IMAGES FOR LM
coneImgEgLminusM = []; % FOR SAVING OUT CONE IMAGES FOR L-M
absorptionsS_LMS = []; % FOR SAVING OUT S-CONE IMAGE
absorptionsLM_LMS = []; % FOR SAVING OUT L+M IMAGE

peakCorrAll = []; % ARRAY FOR STORING SIGNAL QUALITY
for i = 1:length(stimInd) % LOOP OVER STIMULUS INDICES
    for j = 1:length(wv2vis) % LOOP OVER WAVELENGTHS TO VISUALIZE
        % SAVE OUT IMAGES FOR BLUE-YELLOW
        [~, coneImgEg, ~, wave, ~, absorptions] = ...
         ARCwvInFocusConesMeanZsandbox(subjNumAll(subjNumInd),stimInd(i),LMSweights,find(wave==wv2vis(j)),dataPath);
         coneImgEgLMS(:,:,i,j) = coneImgEg;
         absorptionsS_LMS(:,:,i,j) = absorptions(:,:,3);
         absorptionsLM_LMS(:,:,i,j) = absorptions(:,:,1)+absorptions(:,:,2);
        
        % SAVE OUT IMAGES FOR L+M
        [~, coneImgEg, ~, wave, ~, absorptions] = ...
         ARCwvInFocusConesMeanZsandbox(subjNumAll(subjNumInd),stimInd(i),[LMSweights(1) LMSweights(2) 0],find(wave==wv2vis(j)),dataPath); 
        coneImgEgLM(:,:,i,j) = coneImgEg;
        
        % SAVE OUT IMAGES FOR L-M
        [~, coneImgEg, ~, wave, ~, absorptions] = ...
         ARCwvInFocusConesMeanZsandbox(subjNumAll(subjNumInd),stimInd(i),LminusMweights,find(wave==wv2vis(j)),dataPath); 
        coneImgEgLminusM(:,:,i,j) = coneImgEg;        
    end
    % CALCULATE SIGNAL QUALITY VS. WAVELENGTH IN FOCUS FOR BLUE-YELLOW
    [~, ~, ~, ~, peakCorr, ~] = ...
    ARCwvInFocusConesMeanZsandbox(subjNumAll(subjNumInd),stimInd(i),LMSweights,1:101,dataPath);
    peakCorrAll(:,end+1) = peakCorr;
    % CALCULATE SIGNAL QUALITY VS. WAVELENGTH IN FOCUS FOR LUMINANCE
    [~, ~, ~, ~, peakCorr, ~] = ...
    ARCwvInFocusConesMeanZsandbox(subjNumAll(subjNumInd),stimInd(i),LMweights,1:101,dataPath);
    peakCorrAll(:,end+1) = peakCorr;
    % CALCULATE SIGNAL QUALITY VS. WAVELENGTH IN FOCUS FOR RED-GREEN
    [~, ~, ~, ~, peakCorr, ~] = ...
    ARCwvInFocusConesMeanZsandbox(subjNumAll(subjNumInd),stimInd(i),LminusMweights,1:101,dataPath);
    peakCorrAll(:,end+1) = peakCorr;
end
%% PLOT SIGNAL QUALITY

figure; 
plot(wave,peakCorrAll(:,[1 4 7]),'LineWidth',1); 
axis square; 
xlim([400 700]); 
ylim([0.4 1]);
formatFigure('Wavelength (nm)','Signal quality');
legend('1/4','1/1','4/1');
figure; 
plot(wave,peakCorrAll(:,[2 5 8]),'LineWidth',1);
axis square; 
xlim([400 700]); 
ylim([0.4 1]);
formatFigure('Wavelength (nm)','Signal quality');
legend('1/4','1/1','4/1');

%% PLOT DIFFERENT CONE IMAGES FOR BUILDING MODEL INTUITIONS

% SPECIFY WHICH IMAGES TO PLOT
absorptionsS_LMSplot = absorptionsS_LMS;
absorptionsLM_LMSplot = absorptionsLM_LMS;
coneImgEgLMSplot = coneImgEgLMS;
redblue = make_yellow_blue_colormap(1);
axisLims = [min([coneImgEgLM(:); coneImgEgLMS(:)]) ...
            max([coneImgEgLM(:); coneImgEgLMS(:)])];
stimCell = {'More blue' 'Purple' 'More red'};
bColorBar = false;

figure;
set(gcf,'Position',[29 90 890 908]);
for i = 1:length(wv2vis)
    subplot(3,3,(i-1)*3+3);
    imagesc(squeeze(coneImgEgLMSplot(:,:,1,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'LineWidth',2);
    colormap(redblue); 
    xlim([50 120]+57); 
    ylim([55 125]); 
    if i==1 && bColorBar
       colorbar;
    end
    clim(axisLims(2).*[-1 1])
    if bColorBar
        title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{1} ', LMS']);
    end
    subplot(3,3,(i-1)*3+2);
    imagesc(squeeze(absorptionsS_LMSplot(:,:,1,i)));
    axis square; 
    colormap(redblue); 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]); 
    set(gca,'LineWidth',2);
    xlim([50 120]+57); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);
    if bColorBar
       title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{1} ', S image']);
    end
    subplot(3,3,(i-1)*3+1);
    imagesc(squeeze(absorptionsLM_LMSplot(:,:,1,i)));
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]); 
    set(gca,'LineWidth',2);
    colormap(redblue); 
    xlim([50 120]+57); 
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
    imagesc(squeeze(coneImgEgLMSplot(:,:,3,i))); 
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    set(gca,'LineWidth',2);
    colormap(redblue); 
    xlim([50 120]+57); 
    ylim([55 125]); 
    if i==1 && bColorBar
       colorbar;
    end
    clim(axisLims(2).*[-1 1])
    if bColorBar
        title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{3} ', LMS']);
    end
    subplot(3,3,(i-1)*3+2);
    imagesc(squeeze(absorptionsS_LMSplot(:,:,3,i)));
    axis square; 
    colormap(redblue); 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]); 
    set(gca,'LineWidth',2);
    xlim([50 120]+57); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);  
    if bColorBar
        title(['\lambda_{focus} = ' num2str(wv2vis(i)) ', ' stimCell{3} ', S image']);
    end
    subplot(3,3,(i-1)*3+1);
    imagesc(squeeze(absorptionsLM_LMSplot(:,:,3,i)));
    axis square; 
    set(gca,'XTick',[]);
    set(gca,'YTick',[]); 
    set(gca,'LineWidth',2);
    colormap(redblue); 
    xlim([50 120]+57); 
    ylim([55 125]); 
    % colorbar; 
    clim(axisLims(2).*[-1 1]);
    if bColorBar
       title('L+M image');
    end
end