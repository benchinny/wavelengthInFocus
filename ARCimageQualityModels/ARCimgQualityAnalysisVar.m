%% Initialize and clear
ieInit;

%% Set up display struct and build Ben's stimulus

subjNum = 2;

% Setting up display properties
d = displayCreate('OLED-Samsung');
d = displaySet(d, 'name', 'my display');
d = displaySet(d,'ViewingDistance',1); % simulated screen distance

bUseBVAMScal = 1; % if using BVAMS calibration data

if bUseBVAMScal
    drivePath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/BVAMS_calibration_files/display calibration on August3/';
    load([drivePath 'Right_disp_Red.mat']);
    d.spd(:,1) = CurrentSpectrum.Spectral.emission_data;
    load([drivePath 'Right_disp_Green.mat']);
    d.spd(:,2) = CurrentSpectrum.Spectral.emission_data;
    load([drivePath 'Right_disp_Blue.mat']);
    d.spd(:,3) = CurrentSpectrum.Spectral.emission_data;
end
d.gamma(:,1) = (d.gamma(:,1).^(1/2.2)).^2.4;
d.gamma(:,2) = (d.gamma(:,2).^(1/2.2)).^2.6;
d.gamma(:,3) = (d.gamma(:,3).^(1/2.2)).^2.2;

% COLOR MATCHING FUNCTIONS
S = [380 4 101]; % weird convention used by Brainard lab for defining wavelengths
load T_xyz1931; % load color matching functions
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale
wave = S(1):S(2):S(1)+S(2)*(S(3)-1); % define wavelength vector
% DEFOCUSES TO LOOK AT
Dall = -humanWaveDefocus(wave);
% WAVELENGTHS TO LOOK AT
% wvAll = humanWaveDefocusInvert(-Dall);

% PARAMETERS OF WAVEFRONT ANALYSIS
PARAMS.PixelDimension = 512;% size of pupil aperture field in pixels (this defines the resolution of the calculation)
PARAMS.PupilSize = 7; %default values - will be replaced depending on choices below
PARAMS.PupilFieldSize =42; %default values - will be replaced depending on choices below
PARAMS.PupilFitSize = 7; %default values - will be replaced depending on choices below
PARAMS.ImagingWavelength = 0.55;% imaging wavelength in microns
PARAMS.WavefrontResolution = 53;% increase to enhance the display of the wavefront (doesn't affect calculation)

%%

if subjNum==1
    subjName = 'BenChin-OS';
    blockNums = [2 3 4 5 6];
    trialNums = [[1:20]' [1:20]' [1:20]' [1:20]' [1:20]'];
    % blockNums = [2 3];
    % trialNums = [[1:20]' [1:20]']; 
elseif subjNum==2
    subjName = 'S2-OS';
    blockNums = [2 3 4 5 6];
    trialNums = [[1:20]' [1:20]' [1:20]' [1:20]' [1:20]'];
    % blockNums = [2 3];
    % trialNums = [[1:20]' [1:20]'];     
end

wvInFocus1all = zeros([20 1]);
meanv00all = zeros([20 1]);
rgb1all = zeros([20 3]);
defocusBasic = zeros([20 1]);

for l = 5 % LOOP OVER BLOCK
    parfor k = 1:20 % LOOP OVER TRIAL
        % LOADING DATA
        blockNumInd = l;
        blockNumTmp = blockNums(blockNumInd);
        trialNumTmp = trialNums(k,blockNumInd);
        
        AFCp = ARCloadFileBVAMS(subjNum,blockNumTmp); % LOAD BVAMS DATA
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumTmp,0);
        % GET THE TIMESTAMP CORRESPONDING TO THE HALFWAY POINT
        t = seconds(TimeStamp)-min(seconds(TimeStamp));
        tHalfway = max(t)/2;
        tDiffFromHalfway = abs(t-tHalfway);
        [~,indMinT] = min(tDiffFromHalfway);
        FrameStart = (indMinT-29):indMinT; % analyze 30 frames

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(30,65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS = struct;
        PARAMS.PupilSize=mean(table2array(ZernikeTable(FrameStart,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(FrameStart,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        c(:,3:NumCoeffs)=table2array(ZernikeTable(FrameStart,11:width(ZernikeTable)));
        meanC = mean(c,1); % TAKE MEAN OF COEFFICIENTS
        
        % STORE COLORS FOR FIRST AND SECOND STIMULI
        rgb00 = [];
        rgb00(1,:) = AFCp.rgb100(trialNumTmp,:);
        rgb00(2,:) = AFCp.rgb200(trialNumTmp,:);
        % rgb00(:,2) = [0; 0];

        % recreate stimulus
        nDotsI = 320;
        rVal = rgb00(1,1);
        gVal = rgb00(1,2);
        bVal = rgb00(1,3);
        im = AFCwordStimImproved('sea',nDotsI.*[1 1],'green');
        imPatternTmp = squeeze(im(:,:,2));
        imPatternTmp = circshift(imPatternTmp,-15,1);
        I(:,:,3) = bVal.*imresize(imPatternTmp,nDotsI.*[1 1],'nearest');
        I(:,:,2) = gVal.*imresize(imPatternTmp,nDotsI.*[1 1],'nearest');
        I(:,:,1) = rVal.*imresize(imPatternTmp,nDotsI.*[1 1],'nearest');
        I = I./255;
        
        % Turn image into 'scene'
        s = sceneFromFile(I, 'rgb', [], d);  % The display is included here
        % I think this copies the struct into an object
        vcAddObject(s); 
        
        % Turning original stimulus into luminance image
        downScale = 1;
        photonsImgXWorig = RGB2XWFormat(s.data.photons);
        energyImgXWorig = Quanta2Energy(wave',photonsImgXWorig);
        energyImgOrig = XW2RGBFormat(energyImgXWorig,320,320);
        
        lumImgOrig = zeros(size(s.data.photons,1),size(s.data.photons,2));
        for j = 1:length(wave)
            lumImgOrig = lumImgOrig+energyImgOrig(:,:,j).*T_sensorXYZ(2,j).*downScale;
        end
        
        figure; 
        set(gcf,'Position',[289 428 1056 420]);
        subplot(1,3,1);
        plot(d.wave,d.spd(:,1),'r','LineWidth',1.5); hold on;
        plot(d.wave,d.spd(:,2),'g','LineWidth',1.5);
        plot(d.wave,d.spd(:,3),'b','LineWidth',1.5);
        axis square;
        formatFigure('Wavelength (\lambda)','Radiance');
        subplot(1,3,2);
        imagesc(I);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        axis square;
        set(gca,'FontSize',15);
        title('Original');
        subplot(1,3,3);
        plot(s.spectrum.wave,squeeze(s.data.photons(160,160,:)),'-k','LineWidth',1);
        formatFigure('Wavelength (\lambda)','Photons');
        axis square;
        
        %% Computing peak correlation for different wavelengths in focus
        
        peakCorr = [];
        peakPSF = [];
        peakImg = [];
        Dall2 = -humanWaveDefocus(wave(16:66));
        wave2 = wave(16:66);

        for i = 1:length(Dall2)
            zCoeffs = [0 meanC(1:end-1)];
            wvfP = wvfCreate('calc wavelengths', wave, ...
                'measured wavelength', humanWaveDefocusInvert(-Dall2(i)), ...
                'zcoeffs', zCoeffs, 'measured pupil', PARAMS.PupilSize, ...
                'name', sprintf('human-%d', PARAMS.PupilSize),'spatial samples',320);
            wvfP.calcpupilMM = PARAMS.PupilSize;
            wvfP.refSizeOfFieldMM = 42;
            wvfP = wvfSet(wvfP, 'zcoeff', 0, 'defocus');
            
            % Convert to siData format as well as wavefront object
            [siPSFData, wvfP] = wvf2SiPsf(wvfP,'showBar',false,'nPSFSamples',320,'umPerSample',1.5212); 
            oi = wvf2oi(wvfP); % CONVERT TO OPTICS OBJECT
            % oi = oiCreateARC('human',wave,Dall(i)); % create optics
            oi = oiCompute(oi, s); % compute optical image of stimulus
        
            photonsImgXW = RGB2XWFormat(oi.data.photons); % FORMATTING
            energyImgXW = Quanta2Energy(wave,photonsImgXW);
            
            lumImgXW = sum(downScale*bsxfun(@times,energyImgXW,T_sensorXYZ(2,:)),2);
            lumImgXYW = XW2RGBFormat(lumImgXW,400,400);
            lumImg = lumImgXYW(41:360,41:360);

            % COMPUTE MAX CORRELATION
            peakCorr(i) = max(max(xcorr2(lumImgOrig,lumImg)));
            % if ismember(round(humanWaveDefocusInvert(-Dall(i))),[460 520 620])
            %     figure;
            %     set(gcf,'Position',[326 418 924 420]);      
            %     subplot(1,2,1);
            %     imagesc(lumImg); axis square; colormap gray;
            %     subplot(1,2,2);
            %     imagesc(lumImgOrig); axis square; colormap gray;    
            %     title(['wavelength in focus: ' num2str(round(humanWaveDefocusInvert(-Dall(i)))) 'nm, ' ...
            %            'max(xcorr) = ' num2str(peakCorr(i))]);
            % end
            display(['Peak correlation loop ' num2str(i)]);
        end
        
        %% Plotting peak correlation with wavelength in focus
        
        [~,indPeak2] = max(peakCorr);
        wvInFocus1 = wave2(indPeak2);
        wvInFocus1all(k,:) = wvInFocus1;
        rgb1all(k,:) = rgb00(1,:);
        meanv00all(k,:) = AFCp.meanv00(trialNumTmp);
        defocusBasic(k,:) = meanC(4);
        
        figure; 
        hold on;
        % plot(humanWaveDefocusInvert(-Dall),peakCorr./max(peakCorr),'k-','LineWidth',1);
        plot(humanWaveDefocusInvert(-Dall2),peakCorr./max(peakCorr),'k-','LineWidth',1);
        % plot(humanWaveDefocusInvert(-Dall),peakPSF./max(peakPSF),'k-','LineWidth',1);
        % plot(humanWaveDefocusInvert(-Dall),peakImg./max(peakImg),'k-','LineWidth',1);
        axis square;
        set(gca,'FontSize',15);
        xlabel('Wavelength in focus');
        ylabel('Peak correlation');
    end
end

save('/Users/benjaminchin/Documents/ARchromaScraps/ARCmodelOutput2_5.mat','meanv00all','wvInFocus1all','rgb1all','defocusBasic');
