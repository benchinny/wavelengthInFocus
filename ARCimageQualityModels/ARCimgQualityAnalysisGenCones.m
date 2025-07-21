function [wvInFocusST, wvInFocusXC] = ARCimgQualityAnalysisGenCones(stimNum)

%% Initialize and clear
ieInit;

% COLOR CONDITIONS
rgbConditions = [0.555 0.000 1.00; ...
                 0.416 0.320 1.00; ...
                 0.312 0.320 1.00; ...
                 0.555 0.320 0.73; ...
                 0.555 0.320 0.533; ...
                 0.555 0.000 1.00; ...
                 0.416 0.000 1.00; ...
                 0.312 0.000 1.00; ...
                 0.555 0.000 0.73; ...
                 0.555 0.000 0.533; ...
                 ];

%% Set up display struct and build Ben's stimulus

% Setting up display properties
d = displayCreate('OLED-Samsung');
d = displaySet(d, 'name', 'my display');
d = displaySet(d,'ViewingDistance',1); % simulated screen distance
d = displaySet(d,'dpi',378); % simulated screen distance

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

% DEFINING COLOR-MATCHING FUNCTIONS
S = [380 4 101]; % weird convention used by Brainard lab for defining wavelengths
load T_xyz1931; % load color matching functions
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale
wave = S(1):S(2):S(1)+S(2)*(S(3)-1); % define wavelength vector

subjNum = 1;
blockNumTmp = 2;
if subjNum==1
    subjName = 'BenChin-OS';
elseif subjNum==2
    subjName = 'S2-OS'; 
end
trialNumTmp = 1;

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

wvInFocusXC = zeros([size(rgbConditions,1) 1]);
wvInFocusST = zeros([size(rgbConditions,1) 1]);

for k = 1:size(rgbConditions,1)
    % Ben's stimulus
    rVal = rgbConditions(k,1);
    gVal = rgbConditions(k,2);
    bVal = rgbConditions(k,3);

    im1 = imread(['/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/word_image_0' num2str(stimNum) '.png']);
    im1(im1>0) = 255;
    % im1 = flipud(im1);   
    imPatternTmp = squeeze(im1(:,:,3));
    imPatternTmp = [zeros([30 size(imPatternTmp,2)]); imPatternTmp; zeros([30 size(imPatternTmp,2)])];
    imPatternTmp = [zeros([size(imPatternTmp,1) 30]) imPatternTmp zeros([size(imPatternTmp,1) 30])];
    indTestMax = find(imPatternTmp==255);
    [rowTest,colTest]=ind2sub(size(imPatternTmp),indTestMax(1));

    I(:,:,3) = bVal.*double(imPatternTmp);
    I(:,:,2) = gVal.*double(imPatternTmp);
    I(:,:,1) = rVal.*double(imPatternTmp);
    I = I./255;
    
    % Turn image into 'scene'
    s = sceneFromFile(I, 'rgb', [], d);  % The display is included here
    % I think this copies the struct into an object
    vcAddObject(s); 
    % s.data.photons(160,160,:) = ones(size(s.data.photons(160,160,:))).*4e14;
    
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
    plot(s.spectrum.wave,squeeze(s.data.photons(rowTest,colTest,:)),'-k','LineWidth',1);
    formatFigure('Wavelength (\lambda)','Photons');
    axis square;
    
    %% Turning original stimulus into luminance image
    
    downScale = 1;
    photonsImgXWorig = RGB2XWFormat(s.data.photons);
    energyImgXWorig = Quanta2Energy(wave',photonsImgXWorig);
    energyImgOrig = XW2RGBFormat(energyImgXWorig,size(s.data.photons,1),size(s.data.photons,2));
    
    lumImgOrig = zeros(size(s.data.photons,1),size(s.data.photons,2));
    for j = 1:length(wave)
        lumImgOrig = lumImgOrig+energyImgOrig(:,:,j).*T_sensorXYZ(2,j).*downScale;
    end

    %% Computing visual Strehl ratio
    
    peakPSF = [];
    maxRawPSFcheck = [];
    Dall2 = -humanWaveDefocus(wave(16:101));
    peakCorr = [];
    
    for i = 1:length(Dall2)
        zCoeffs = [0 meanC(1:end-1)];
        wvfP = wvfCreate('calc wavelengths', wave, ...
            'measured wavelength', humanWaveDefocusInvert(-Dall2(i)), ...
            'zcoeffs', zCoeffs, 'measured pupil', PARAMS.PupilSize, ...
            'name', sprintf('human-%d', PARAMS.PupilSize),'spatial samples',max(size(imPatternTmp)));
        wvfP.calcpupilMM = PARAMS.PupilSize;
        wvfP.refSizeOfFieldMM = 42;
        wvfP = wvfSet(wvfP, 'zcoeff', 0, 'defocus');
        
        % Convert to siData format as well as wavefront object
        [siPSFData, wvfP] = wvf2SiPsf(wvfP,'showBar',false,'nPSFSamples',max(size(imPatternTmp)),'umPerSample',1.1512); 
        oi = wvf2oi(wvfP); % CONVERT TO OPTICS OBJECT
        % oi = oiCreateARC('human',wave,Dall(i)); % create optics
        oi = oiCompute(oi, s); % compute optical image of stimulus

        % Create the coneMosaic object
        cMosaic = coneMosaic;

        % Set size to show relevant portion of scene
        cMosaic.setSizeToFOV(1 * sceneGet(s, 'fov'));

        % key line for computing absorptions
        absorptions = cMosaic.computeSingleFrame(oi, 'fullLMS', true); 
        
        maxRawPSFcheck(i) = max(max(ifftshift(ifft2(oi.optics.OTF.OTF(:,:,i)))));
    
        photonsPSFXW = RGB2XWFormat(siPSFData.psf); % FORMATTING
        photonsPSFXWweighted = bsxfun(@times,photonsPSFXW,squeeze(s.data.photons(rowTest,colTest,:))');
        energyPSFXW = Quanta2Energy(wave,photonsPSFXWweighted);
        
        lumPSFXW = sum(bsxfun(@times,energyPSFXW,T_sensorXYZ(2,:)),2);
        peakPSF(i) = max(lumPSFXW);

        photonsImgXW = RGB2XWFormat(oi.data.photons);
        energyImgXW = Quanta2Energy(wave,photonsImgXW);
        energyImg = XW2RGBFormat(energyImgXW,size(oi.data.photons,1),size(oi.data.photons,2));
        
        lumImg = zeros(size(oi.data.photons,1),size(oi.data.photons,2));
        for j = 1:length(wave)
           lumImg = lumImg+downScale*energyImg(:,:,j).*T_sensorXYZ(2,j);
        end
        cropCorner = floor((size(lumImg)-size(lumImgOrig))/2);
        lumImg = imcrop(lumImg,[fliplr(cropCorner) fliplr(size(lumImgOrig))]);
        peakCorr(i) = max(max(xcorr2(lumImgOrig,lumImg)));
        % if ismember(round(humanWaveDefocusInvert(-Dall2(i))),[460 520 620])
        %     figure;
        %     set(gcf,'Position',[326 418 924 420]);      
        %     subplot(1,2,1);
        %     imagesc(lumImg); axis square; colormap gray;
        %     subplot(1,2,2);
        %     imagesc(lumImgOrig); axis square; colormap gray;    
        %     title(['wavelength in focus: ' num2str(round(humanWaveDefocusInvert(-Dall2(i)))) 'nm, ' ...
        %            'max(xcorr) = ' num2str(peakCorr(i))]);
        % end
        display(['Correlation iteration ' num2str(i)]);        
    end
end

