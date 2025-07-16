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
nDotsI = 260;
rVal = 0.00;
bVal = 1.00;
im = AFCwordStimImproved('sea',nDotsI.*[1 1],'green');
imPatternTmp = squeeze(im(:,:,2));
imPatternTmp = circshift(imPatternTmp,-15,1);
I(:,:,3) = bVal.*imresize(imPatternTmp,nDotsI.*[1 1],'nearest');
I(:,:,1) = rVal.*imresize(imPatternTmp,nDotsI.*[1 1],'nearest');
I = I./255;
I = zeros(size(I));
% I(159:161,159:161,1) = rVal;
% I(159:161,159:161,2) = 0.00;
% I(159:161,159:161,3) = bVal;
I(130,130,1) = rVal;
I(130,130,2) = 0.00;
I(130,130,3) = bVal;
% I(156:164,156:164,1) = rVal;
% I(156:164,156:164,2) = 0.00;
% I(156:164,156:164,3) = bVal;
% I(147:173,147:173,1) = rVal;
% I(147:173,147:173,2) = 0.00;
% I(147:173,147:173,3) = bVal;

% GABOR
frqCpd = 15;
contrast = 1;
rgbAll = [rVal 0 bVal];
k0 = 1;
gammaR = 2.4;
gammaG = 2.6;
gammaB = 2.2;

acuStimOrig1 = ARC2Dgabor(smpPos(260,260),[],0,0,[frqCpd 3*frqCpd 5*frqCpd 7*frqCpd], ...
               [contrast contrast/3 contrast/5 contrast/7],-15,90,0.2,0.2, ...
               [rgbAll(k0,1)^gammaR rgbAll(k0,2)^gammaG rgbAll(k0,3)^gammaB],1,1,0,0);

acuStimOrig1(:,:,1) = acuStimOrig1(:,:,1).^(1/gammaR);
acuStimOrig1(:,:,2) = acuStimOrig1(:,:,2).^(1/gammaG);
acuStimOrig1(:,:,3) = acuStimOrig1(:,:,3).^(1/gammaB);
I = acuStimOrig1;

% Turn image into 'scene'
s = sceneFromFile(I, 'rgb', [], d);  % The display is included here
% I think this copies the struct into an object
vcAddObject(s); 
% s.data.photons(130,130,:) = ones(size(s.data.photons(130,130,:))).*4e14;

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
plot(s.spectrum.wave,squeeze(s.data.photons(130,130,:)),'-k','LineWidth',1);
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
    
    photonsTmp = squeeze(s.data.photons(130,130,:));
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
energyImgOrig = XW2RGBFormat(energyImgXWorig,260,260);

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
    energyImg = XW2RGBFormat(energyImgXW,326,326);
    
    lumImg = zeros(size(oi.data.photons,1),size(oi.data.photons,2));
    for j = 1:length(wave)
       lumImg = lumImg+downScale*energyImg(:,:,j).*T_sensorXYZ(2,j);
    end
    lumImg = lumImg(33:292,33:292);
%    lumImg = oi.data.illuminance(41:360,41:360);
    % lumImg = lumImg(33:288,33:288);
    peakCorr(i) = max(max(xcorr2(lumImgOrig,lumImg)));
    if ismember(round(humanWaveDefocusInvert(-Dall2(i))),[435 480 550])
        figure;
        set(gcf,'Position',[326 418 924 420]);      
        subplot(1,2,1);
        imagesc(lumImg); axis square; colormap gray;
        subplot(1,2,2);
        imagesc(lumImgOrig); axis square; colormap gray;    
        title(['wavelength in focus: ' num2str(round(humanWaveDefocusInvert(-Dall2(i)))) 'nm, ' ...
               'max(xcorr) = ' num2str(peakCorr(i))]);
    end
    display(['Iteration ' num2str(i)]);
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
