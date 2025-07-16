%% SUBJECTS AND STIMULI TO ANALYZE

subjNumAll = [1 3 5 10 16 17 18 20];
stimNumAll = [1 3 5];
stimNumReal = [1 8 6];
imgQualQueryAll = [];
wvInFocusRawAll = [];
imgQualAll = [];
wLMSall = [];
lagLeadAll = [];

%% GET IMAGE QUALITY VALUES

for i = 1:length(stimNumAll)
    for j = 1:length(subjNumAll)
        [waveInFocus3raw, waveInFocus3, wLMS, lagLead] = ...
        ARCtestWvInFocusMeanZspatFilterLMSlagLeadEffect(subjNumAll(j),'LMS',stimNumAll(i));
        [imgQualQuery, peakCorr] = ARCwvInFocusConesMeanZspatFilterImgQual(subjNumAll(j),stimNumReal(i),wLMS,[waveInFocus3raw waveInFocus3]);
        imgQualQueryAll(i,j,:) = imgQualQuery;
        wvInFocusRawAll(i,j,:) = [waveInFocus3raw waveInFocus3];
        imgQualAll(i,j,:) = peakCorr;
        wLMSall(j,:) = wLMS;
        lagLeadAll(i,j,:) = lagLead;
    end
end

%% PLOTTING

load('/Users/benjaminchin/Documents/ARchromaScraps/imgQualAll.mat');
load('/Users/benjaminchin/Documents/ARchromaScraps/imgQualQueryAll.mat');
load('/Users/benjaminchin/Documents/ARchromaScraps/wvInFocusRawAll.mat');
load('/Users/benjaminchin/Documents/ARchromaScraps/wLMSall.mat');

stimString = {'More blue' 'Purple' 'More red'};
wave = 380:4:780;

for j = 1:length(stimNumAll)
    figure;
    set(gcf,'Position',[311 330 1236 594]);
    for i = 1:length(subjNumAll)
        subplot(2,4,i);
        hold on;
        plot(wave,squeeze(imgQualAll(j,i,:)),'k-','LineWidth',1);
        plot(squeeze(wvInFocusRawAll(j,i,2:4)),squeeze(imgQualQueryAll(j,i,2:4)),'kx-','LineWidth',1,'MarkerSize',10);
        plot(squeeze(wvInFocusRawAll(j,i,1)),squeeze(imgQualQueryAll(j,i,1)),'ko','LineWidth',1,'MarkerFaceColor','w','MarkerSize',10);
        axis square;
        ylim([0 1]);
        xlim([400 700]);
        set(gca,'FontSize',15);
        ylabel('Peak cross-correlation');
        xlabel('Wavelength in focus');
        if i==1
            title(stimString{j});
        end
    end
end

%% VISUALIZATION FOR CONE MECHANISM

subjInd = [1 6];
stimString = {'More blue' 'Purple' 'More red'};

for i = 1:length(subjInd)
    figure;
    set(gcf,'Position',[235 381 1241 558]);
    for j = 1:size(wvInFocusRawAll,1)
        imgQuery = ARCwvInFocusConesMeanZspatFilterImgQuery(subjNumAll(subjInd(i)),stimNumReal(j),wLMSall(subjInd(i),:),squeeze(wvInFocusRawAll(j,subjInd(i),:)));
        wvInFocusTmp = flipud(squeeze(wvInFocusRawAll(j,subjInd(i),:)));
        for k = 1:size(wvInFocusRawAll,3)
            subplot(3,4,(j-1)*4+k);
            imagesc(squeeze(imgQuery(:,:,k)));
            axis square;
            colormap gray;
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            title(['Wavelength = ' num2str(wvInFocusTmp(k),3) ', ' stimString{j}]);
        end
    end
end

%% VISUALIZATION FOR CONE MECHANISM

subjInd = [1 6];
stimString = {'More blue' 'Purple' 'More red'};

for i = 1:length(subjInd)
    figure;
    set(gcf,'Position',[235 381 1241 558]);
    for j = 1:size(wvInFocusRawAll,1)
        imgQuery = ARCwvInFocusConesMeanZspatFilterImgQuery(subjNumAll(subjInd(i)),stimNumReal(j),[0.72 0.28 0],squeeze(wvInFocusRawAll(j,subjInd(i),:)));
        wvInFocusTmp = flipud(squeeze(wvInFocusRawAll(j,subjInd(i),:)));
        for k = 1:size(wvInFocusRawAll,3)
            subplot(3,4,(j-1)*4+k);
            imagesc(squeeze(imgQuery(:,:,k)));
            axis square;
            colormap gray;
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            title(['Wavelength = ' num2str(wvInFocusTmp(k),3) ', ' stimString{j}]);
        end
    end
end


%%

defocusAtRef = -[0:0.1:4];
penalty = defocusAtRef.^2;
stimDist = 3.5;
wvInFocusPossible = 400:700;
wvRef = 550;

figure;
plot(defocusAtRef,penalty,'k-','LineWidth',1)
axis square;
set(gca,'FontSize',18);
xlabel('Defocus at 550nm');
ylabel('Penalty');

figure; 
plot(wvInFocusPossible,(humanWaveDefocusARC(wvInFocusPossible,wvRef,1)-stimDist).^2,'k-','LineWidth',1);
axis square;
set(gca,'FontSize',18);
xlabel('Wavelength-in-focus (nm)');
ylabel('Penalty (unitless?)');

%%

estAcc = [];

for i = 1:8
    estAcc(:,i) = [1.5 2.5 3.5]'-squeeze(lagLeadAll(1,i,:))+1;
end
