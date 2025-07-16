%% Initialize and clear
ieInit;

% COLOR CONDITIONS
rgbConditions = [0.555 0.320 1.00; ...
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
% d.spd = ones(size(d.spd)).*0.000200 + 0.03.*(repmat([1:length(d.wave)]',[1 3])./1000);
% d.spd(:,1) = [normpdf(380:4:780,624,10).*0.005]';
% d.spd(:,2) = [normpdf(380:4:780,532,10).*0.005]';
% d.spd(:,3) = [normpdf(380:4:780,488,10).*0.005]';

% DEFINING COLOR-MATCHING FUNCTIONS
S = [380 4 101]; % weird convention used by Brainard lab for defining wavelengths
load T_xyz1931; % load color matching functions
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale
wave = S(1):S(2):S(1)+S(2)*(S(3)-1); % define wavelength vector
% T_sensorXYZ(2,:) = normpdf(wave,556,40).*70000;

% MAKING 2D CSF FUNCTION
oi = oiCreateARC('human',wave,0); % create optics

[fx, fy] = meshgrid(oi.optics.OTF.fx,oi.optics.OTF.fy);
% % scale so frequencies are in units of cyc/deg
fx = fx./3.37;
fy = fy./3.37;
df = sqrt(fx.^2 + fy.^2); % compute distance from origin
CSF2d = 0.04992*(1+5.9375*df).*exp(-0.114*df.^1.1);
% inverse Fourier transform of 2D CSF
N = ifftshift(ifft2(fftshift(CSF2d)));

for k = 1:size(rgbConditions,1)
    % Ben's stimulus
    rVal = rgbConditions(k,1);
    gVal = rgbConditions(k,2);
    bVal = rgbConditions(k,3);

    im1 = imread('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/word_image_04.png');
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
    
    %% Computing visual Strehl ratio
    
    Dall = -humanWaveDefocus(wave); % defocus values to look at
    peakPSF = [];
    polyPSFall = [];
    maxRawPSFcheck = [];
    
    for i = 1:length(Dall)
    
        oi = oiCreateARC('human',wave,Dall(i)); % create optics
        maxRawPSFcheck(i) = max(max(ifftshift(ifft2(oi.optics.OTF.OTF(:,:,i)))));
    
        polyPSF = [];
        
        photonsTmp = squeeze(s.data.photons(rowTest,colTest,:));
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
    [~,ind] = max(peakPSF);
    wvInFocusST(k) = humanWaveDefocusInvert(-Dall(ind));
    % original: 612 528 476 612 612 616 472 472 620 620
    % ave:      612 528 476 612 612 616 472 472 620 620
    % ceo:      612 528 476 612 612 616 472 472 620 620
    % sea:      612 528 476 612 612 616 472 472 620 620
    % osx:      612 528 476 612 612 616 472 472 620 620
    %% Turning original stimulus into luminance image
    
    downScale = 1;
    photonsImgXWorig = RGB2XWFormat(s.data.photons);
    energyImgXWorig = Quanta2Energy(wave',photonsImgXWorig);
    energyImgOrig = XW2RGBFormat(energyImgXWorig,size(s.data.photons,1),size(s.data.photons,2));
    
    lumImgOrig = zeros(size(s.data.photons,1),size(s.data.photons,2));
    for j = 1:length(wave)
        lumImgOrig = lumImgOrig+energyImgOrig(:,:,j).*T_sensorXYZ(2,j).*downScale;
    end
    
    %% Computing peak correlation for different wavelengths in focus
    
    peakCorr = [];
    % Dall2 = fliplr(Dall);
    Dall2 = -humanWaveDefocus(wave);
    peakPSF = [];
    peakImg = [];
    
    for i = 1:length(Dall2)
    
        oi = oiCreateARC('human',wave,Dall2(i)); % create optics
        oi = oiCompute(oi, s); % compute optical image of stimulus
    
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
        if ismember(round(humanWaveDefocusInvert(-Dall2(i))),[460 520 620])
            figure;
            set(gcf,'Position',[326 418 924 420]);      
            subplot(1,2,1);
            imagesc(lumImg); axis square; colormap gray;
            subplot(1,2,2);
            imagesc(lumImgOrig); axis square; colormap gray;    
            title(['wavelength in focus: ' num2str(round(humanWaveDefocusInvert(-Dall2(i)))) 'nm, ' ...
                   'max(xcorr) = ' num2str(peakCorr(i))]);
        end
        display(['Correlation iteration ' num2str(i)]);
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

    [~,ind] = max(peakCorr);
    wvInFocusXC(k) = humanWaveDefocusInvert(-Dall2(ind));
    % original: 552 532 520 572 584 564 516 496 604 612
    % ave:      564 532 520 580 588 592 496 488 612 616
    % ceo:      564 532 520 580 588 592 496 488 612 616
    % sea:      564 532 520 580 588 596 496 488 612 616
    % osx:      564 532 520 580 588 592 496 488 612 616
end

