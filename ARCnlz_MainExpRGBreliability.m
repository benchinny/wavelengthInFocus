%% LOOKING AT RELIABILITY OF WAVEFRONT MEASUREMENTS

subjNum = 11;
plot2make = 'raw';
% plot2make = 'correlation';

if subjNum==11 || subjNum==12
   blockNums = 2:7;
   trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   subjName = ['S' num2str(subjNum) '-OD'];
elseif subjNum==13
   blockNums = 3:8;
   trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   subjName = ['S' num2str(subjNum) '-OD'];   
end

rgb1all = [];
meanv00all = [];

for i = 1:length(blockNums)
    blockNumTmp = blockNums(i);
    trialNumsTmp = trialNums{i};
    AFCp = ARCloadFileBVAMS(subjNum,blockNumTmp);
    for j = 1:length(trialNumsTmp)
        [ZernikeTable, ~, ~, ~] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumsTmp(j),0);
        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65.
        
        PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(:,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0;
        if strcmp(plot2make,'raw')
            figure(1); 
            plot(c(:,1:7)); 
            axis square;
            set(gca,'FontSize',15);
            xlabel('Frame');
            ylabel('Coefficient');
            title(['Block ' num2str(blockNumTmp) ' Trial ' num2str(j)]);
        end
        if strcmp(plot2make,'correlation')
            figure(1);
            set(gcf,'Position',[233 291 1363 631]);
            subplot(2,3,1);
            plot(c(~indBad,4),c(~indBad,3),'ko');
            title(['Correlation = ' num2str(corr(c(~indBad,4),c(~indBad,3)),3)]);
            axis square;
            set(gca,'FontSize',15);
            xlabel('Coefficient 4 (defocus)');
            ylabel('Coefficient 3 (oblique astigmatism)');
            subplot(2,3,2);
            plot(c(~indBad,4),c(~indBad,5),'ko');
            title(['Correlation = ' num2str(corr(c(~indBad,4),c(~indBad,5)),3)]);
            axis square;
            set(gca,'FontSize',15);
            xlabel('Coefficient 4 (defocus)');
            ylabel('Coefficient 5 (astigmatism)');       
            subplot(2,3,3);
            plot(c(~indBad,4),c(~indBad,6),'ko');
            title(['Correlation = ' num2str(corr(c(~indBad,4),c(~indBad,6)),3)]);
            axis square;
            set(gca,'FontSize',15);
            xlabel('Coefficient 4 (defocus)');
            ylabel('Coefficient 6 (oblique trefoil)');       
            subplot(2,3,4);
            plot(c(~indBad,4),c(~indBad,7),'ko');
            title(['Correlation = ' num2str(corr(c(~indBad,4),c(~indBad,7)),3)]);
            axis square;
            set(gca,'FontSize',15);
            xlabel('Coefficient 4 (defocus)');
            ylabel('Coefficient 7 (vertical coma)');     
            subplot(2,3,5);
            plot(c(~indBad,4),c(~indBad,8),'ko');
            title(['Correlation = ' num2str(corr(c(~indBad,4),c(~indBad,8)),3)]);
            axis square;
            set(gca,'FontSize',15);
            xlabel('Coefficient 4 (defocus)');
            ylabel('Coefficient 8 (horizontal coma)');      
            subplot(2,3,6);
            plot(c(~indBad,4),c(~indBad,9),'ko');
            title(['Correlation = ' num2str(corr(c(~indBad,4),c(~indBad,9)),3)]);
            axis square;
            set(gca,'FontSize',15);
            xlabel('Coefficient 4 (defocus)');
            ylabel('Coefficient 9 (horizontal trefoil)');                  
        end
        pause;
    end
    rgb1all = [rgb1all; AFCp.rgb100];
    meanv00all = [meanv00all; AFCp.meanv00./1.2255];
end

%% LOOKING AT PSFs AT EVERY TIME POINT

subjNum = 12;
subjName = ['S' num2str(subjNum) '-OD'];
plot2make = 'raw';
% plot2make = 'correlation';

blockNumTmp = 2;
trialNumsTmp = 11;
AFCp = ARCloadFileBVAMS(subjNum,blockNumTmp);
[ZernikeTable, ~, ~, ~] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumsTmp,0);
NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65.

PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
PARAMS.PupilFitSize=mean(table2array(ZernikeTable(:,5))); 
PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
indGd = find(c(:,4)~=0);

wave = 380:5:780;
for j = 1:length(indGd)
    zCoeffs = [0 c(j,1:end-1)];
    wvfP = wvfCreate('calc wavelengths', 875, ...
        'measured wavelength', 875, ...
        'zcoeffs', zCoeffs, 'measured pupil', PARAMS.PupilSize, ...
        'name', sprintf('human-%d', PARAMS.PupilSize),'spatial samples',320);
    wvfP.calcpupilMM = PARAMS.PupilSize;
    wvfP.refSizeOfFieldMM = 42;
    wvfP = wvfSet(wvfP, 'zcoeff', 0, 'defocus');
    
    % Convert to siData format as well as wavefront object
    [siPSFData, wvfP] = wvf2SiPsf(wvfP,'showBar',false,'nPSFSamples',320,'umPerSample',1.1512);
    figure(1);
    imagesc(squeeze(siPSFData.psf(:,:,1)));
    axis square;
    pause;
end
