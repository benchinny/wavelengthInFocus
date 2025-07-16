%% Initialize and clear
ieInit;

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
% d.spd = ones(size(d.spd)).*0.000200 + 0.03.*(repmat([1:length(d.wave)]',[1 3])./1000);
% d.spd(:,1) = [normpdf(380:4:780,624,10).*0.005]';
% d.spd(:,2) = [normpdf(380:4:780,532,10).*0.005]';
% d.spd(:,3) = [normpdf(380:4:780,488,10).*0.005]';

% Ben's stimulus
nDotsI = 320;
rVal = 0.569;
gVal = 0.432;
bVal = 1.00;
im = imread('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/word_image_01.png');
im = double(im);
imPatternTmp = squeeze(im(:,:,3));
I(:,:,3) = bVal.*imresize(imPatternTmp,nDotsI.*[1 1],'nearest');
I(:,:,2) = gVal.*imresize(imPatternTmp,nDotsI.*[1 1],'nearest');
I(:,:,1) = rVal.*imresize(imPatternTmp,nDotsI.*[1 1],'nearest');
I = I./255;
% I(159:161,159:161,1) = rVal;
% I(159:161,159:161,2) = 0.00;
% I(159:161,159:161,3) = bVal;
% I(160,160,1) = rVal;
% I(160,160,2) = 0.00;
% I(160,160,3) = bVal;
% I(156:164,156:164,1) = rVal;
% I(156:164,156:164,2) = 0.00;
% I(156:164,156:164,3) = bVal;
% I(147:173,147:173,1) = rVal;
% I(147:173,147:173,2) = 0.00;
% I(147:173,147:173,3) = bVal;

% acuStimOrig1 = ARC2Dgabor(smpPos(256,256),[],0,0,24,1,-15,0,0.2,0.2,[0.25 0 0],1,1,0,1);
% acuStimOrig1(:,:,1) = acuStimOrig1(:,:,1).^(1/2.4);
% acuStimOrig1(:,:,3) = acuStimOrig1(:,:,3).^(1/2.6);
% I = acuStimOrig1;

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
plot(s.spectrum.wave,squeeze(s.data.photons(160,160,:)),'-k','LineWidth',1);
formatFigure('Wavelength (\lambda)','Photons');
axis square;

%% Computing visual Strehl ratio

S = [380 4 101]; % weird convention used by Brainard lab for defining wavelengths
load T_xyz1931; % load color matching functions
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale
wave = S(1):S(2):S(1)+S(2)*(S(3)-1); % define wavelength vector
% T_sensorXYZ(2,:) = normpdf(wave,556,40).*70000;

oi = oiCreateARC('human',wave,0); % create optics

% % making 2D CSF function
[fx, fy] = meshgrid(oi.optics.OTF.fx,oi.optics.OTF.fy);
% % scale so frequencies are in units of cyc/deg
fx = fx./3.37;
fy = fy./3.37;
df = sqrt(fx.^2 + fy.^2); % compute distance from origin
CSF2d = 0.04992*(1+5.9375*df).*exp(-0.114*df.^1.1);
% inverse Fourier transform of 2D CSF
N = ifftshift(ifft2(fftshift(CSF2d)));

Dall = -humanWaveDefocus(380:5:780); % defocus values to look at
peakPSF = [];
polyPSFall = [];
maxRawPSFcheck = [];

for i = 1:length(Dall)

    oi = oiCreateARC('human',wave,Dall(i)); % create optics
    maxRawPSFcheck(i) = max(max(ifftshift(ifft2(oi.optics.OTF.OTF(:,:,i)))));

    polyPSF = [];
    
    photonsTmp = squeeze(s.data.photons(160,160,:));
    energyTmp = Quanta2Energy(wave,photonsTmp);

    for ind = 1:length(wave)
        polyPSF(:,:,ind) = ifftshift(ifft2(oi.optics.OTF.OTF(:,:,ind))).* ...
                          energyTmp(ind).* ...
                          T_sensorXYZ(2,ind);
    end
    
    polyPSF = sum(polyPSF,3);
    peakPSF(i) = max(max(polyPSF));
    vsx(i) = sum(sum(N.*polyPSF));

    polyPSFall(:,:,i) = polyPSF;
end

figure; 
% plot(humanWaveDefocusInvert(-Dall),vsx./max(vsx),'k-','LineWidth',1); hold on;
plot(humanWaveDefocusInvert(-Dall),0.29.*peakPSF./max(peakPSF),'k-','LineWidth',1);
% legend('Visual Strehl','Strehl');
axis square;
set(gca,'FontSize',15);
xlabel('Wavelength in focus');
ylabel('Ratio');
% ind = 21; % examine at particular wavelength index
% testWave = oi.optics.OTF.wave(ind);
% testOTF = fftshift(oi.optics.OTF.OTF(:,:,ind));
% testPSF = ifftshift(ifft2(oi.optics.OTF.OTF(:,:,ind)));

%% Turning original stimulus into luminance image

downScale = 1;
photonsImgXWorig = RGB2XWFormat(s.data.photons);
energyImgXWorig = Quanta2Energy(wave',photonsImgXWorig);
energyImgOrig = XW2RGBFormat(energyImgXWorig,320,320);

lumImgOrig = zeros(size(s.data.photons,1),size(s.data.photons,2));
for j = 1:length(wave)
    lumImgOrig = lumImgOrig+energyImgOrig(:,:,j).*T_sensorXYZ(2,j).*downScale;
end

%% Computing peak correlation for different wavelengths in focus

peakCorr = [];
% Dall2 = fliplr(Dall);
Dall2 = -humanWaveDefocus(380:5:780);
peakPSF = [];
peakImg = [];

for i = 1:length(Dall2)

    oi = oiCreateARC('human',wave,Dall2(i)); % create optics
    oi = oiCompute(oi, s); % compute optical image of stimulus

    photonsImgXW = RGB2XWFormat(oi.data.photons);
    energyImgXW = Quanta2Energy(wave,photonsImgXW);
    energyImg = XW2RGBFormat(energyImgXW,400,400);
    
    lumImg = zeros(size(oi.data.photons,1),size(oi.data.photons,2));
    for j = 1:length(wave)
       lumImg = lumImg+downScale*energyImg(:,:,j).*T_sensorXYZ(2,j);
    end
    lumImg = lumImg(41:360,41:360);
%    lumImg = oi.data.illuminance(41:360,41:360);
    % lumImg = lumImg(33:288,33:288);
    peakCorr(i) = max(max(xcorr2(lumImgOrig,lumImg)));
    if ismember(round(humanWaveDefocusInvert(-Dall2(i))),[460 520 620])
        if ismember(round(humanWaveDefocusInvert(-Dall2(i))),[620])
            oiIllustrate = oi;
            imgWaveR = zeros([size(oiIllustrate.data.photons,1) size(oiIllustrate.data.photons,2) 3]);
            imgWaveR(:,:,1) = squeeze(oiIllustrate.data.photons(:,:,61))./max(max(squeeze(oiIllustrate.data.photons(:,:,61))));
            imgWaveG = zeros([size(oiIllustrate.data.photons,1) size(oiIllustrate.data.photons,2) 3]);
            imgWaveG(:,:,2) = squeeze(oiIllustrate.data.photons(:,:,38))./max(max(squeeze(oiIllustrate.data.photons(:,:,38))));
            imgWaveGB = zeros([size(oiIllustrate.data.photons,1) size(oiIllustrate.data.photons,2) 3]);
            imgWaveGB(:,:,3) = squeeze(oiIllustrate.data.photons(:,:,30))./max(max(squeeze(oiIllustrate.data.photons(:,:,30))));  
            imgWaveGB(:,:,2) = squeeze(oiIllustrate.data.photons(:,:,30))./max(max(squeeze(oiIllustrate.data.photons(:,:,30))));  
            imgWaveB = zeros([size(oiIllustrate.data.photons,1) size(oiIllustrate.data.photons,2) 3]);
            imgWaveB(:,:,3) = squeeze(oiIllustrate.data.photons(:,:,21))./max(max(squeeze(oiIllustrate.data.photons(:,:,21))));
            imgWaveRG = zeros([size(oiIllustrate.data.photons,1) size(oiIllustrate.data.photons,2) 3]);
            imgWaveRG(:,:,1) = squeeze(oiIllustrate.data.photons(:,:,51))./max(max(squeeze(oiIllustrate.data.photons(:,:,51))));  
            imgWaveRG(:,:,2) = 0.7.*squeeze(oiIllustrate.data.photons(:,:,51))./max(max(squeeze(oiIllustrate.data.photons(:,:,51))));              
            figure;
            set(gcf,'Position',[263 471 1219 420]);
            subplot(3,2,1);
            imagesc(imgWaveR);
            subplot(3,2,2);
            imagesc(imgWaveRG);            
            subplot(3,2,3);
            imagesc(imgWaveG);
            subplot(3,2,4);
            imagesc(imgWaveGB);            
            subplot(3,2,5);
            imagesc(imgWaveB);            
        end
        figure;
        set(gcf,'Position',[326 418 924 420]);      
        subplot(1,2,1);
        imagesc(lumImg); axis square; colormap gray;
        subplot(1,2,2);
        imagesc(lumImgOrig); axis square; colormap gray;    
        title(['wavelength in focus: ' num2str(round(humanWaveDefocusInvert(-Dall2(i)))) 'nm, ' ...
               'max(xcorr) = ' num2str(peakCorr(i))]);
    end
end

%% Plotting peak correlation with wavelength in focus

figure; 
hold on;
% plot(humanWaveDefocusInvert(-Dall2),peakCorr./max(peakCorr),'k-','LineWidth',1);
plot(humanWaveDefocusInvert(-Dall2),peakCorr./max(peakCorr),'k-','LineWidth',1);
% plot(humanWaveDefocusInvert(-Dall2),peakPSF./max(peakPSF),'k-','LineWidth',1);
% plot(humanWaveDefocusInvert(-Dall2),peakImg./max(peakImg),'k-','LineWidth',1);
axis square;
set(gca,'FontSize',15);
xlabel('Wavelength in focus');
ylabel('Peak correlation');

%% MAKING 1D CSF EQUATION

f = oi.optics.OTF.fx(31:60)./3.37; % defining frequency space
CSF1d = 0.04992*(1+5.9375*f).*exp(-0.114*f.^1.1); % 1D CSF equation

%% MAKING 1D CSF EQUATION

[fx, fy] = meshgrid(smpFrq(60,202),smpFrq(60,202)); % defining frequency space
df = sqrt(fx.^2 + fy.^2); % compute distance from origin
CSF2d = 0.04992*(1+5.9375*df).*exp(-0.114*df.^1.1);
% inverse Fourier transform of 2D CSF
N = ifftshift(ifft2(fftshift(CSF2d)));

%% GET DIFFRACTION-LIMITED PSF

oi = oiCreate('diffraction limited'); 
oi = oiSet(oi,'fnumber',5.6); 
% oiPlot(oi,'psf 550');

if notDefined('oi'), oi = vcGetObject('oi'); end
if notDefined('pType'), pType = 'otf550'; end

wavelength = oiGet(oi, 'wavelength');
optics = oiGet(oi, 'optics');

% This catches the case in which the oi has not yet been defined, but the
% optics have.
if isempty(wavelength)
    oi = initDefaultSpectrum(oi, 'hyperspectral');
    optics = initDefaultSpectrum(optics, 'hyperspectral');
    wavelength = opticsGet(optics, 'wavelength');
end

nWave = oiGet(oi, 'nwave');
pType = ieParamFormat(pType);

% Spatial scale is microns.
units = 'um';
nSamp = 100;
freqOverSample = 4;
if strfind(pType, '550')
    thisWave = 550;
elseif length(varargin) >= 1
    thisWave = varargin{1};
else
    thisWave = ieReadNumber('Select PSF wavelength (nm)', 550, ...
        '%.0f');
end

opticsModel = opticsGet(optics, 'model');

thisWave = 400;
psf = opticsGet(optics, 'psf data', ...
                    thisWave, units, nSamp, freqOverSample);

%%

lumRed = [0.25 0.5 0.75 1];
lumBlue = [0.25 0.5 0.75 1];
bestFocusLambdaFullRedCorr = [600 600 580 570];
bestFocusLambdaFullBlueCorr = [500 520 545 570];

figure;
hold on;
plot(lumRed,bestFocusLambdaFullBlueCorr,'-ro','LineWidth',1.5,'MarkerSize',12);
plot(lumBlue,bestFocusLambdaFullRedCorr,'-bo','LineWidth',1.5,'MarkerSize',12);
axis square;
formatFigure('Relative Luminance','Wavelength in focus','Correlation metric');
ylim([400 700]);


