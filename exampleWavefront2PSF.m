%% SCRIPT FOR PLAYING AROUND WITH POINT-SPREAD FUNCTIONS

wave = 555; % WAVELENGTHS TO GENERATE POINT-SPREAD FUNCTION AT
% MAKE SOME EXAMPLE ZERNIKE COEFFICIENTS
zCoeffs = zeros([1 65]);
% NOTE THAT THE FIRST 3 ZERNIKE POLYNOMIALS DO NOT AFFECT THE PSF SHAPE. 4
% 5 AND 6 COVER ASTIGMATISM AND DEFOCUS
zCoeffs(4) = 0.1; % OBLIQUE ASTIGMATISM
zCoeffs(5) = -0.4; % DEFOCUS TERM
zCoeffs(6) = 0.1; % VERTICAL ASTIGMATISM
zCoeffs(7:11) = 0.1.*rand([1 5]);

pupilSize = 4; % PUPIL SIZE
refSizeOfFieldMM = 11; % SAMPLING OF WAVEFRONT (MESS AROUND WITH THIS NUMBER IF PSF LOOKS FUNNY)
nSamples = 260; % NUMBER OF SAMPLES OF BOTH WAVEFRONT AND POINT-SPREAD FUNCTION
wvfP = wvfCreate('calc wavelengths', wave, ...
    'measured wavelength', wave, ... % MAKE SURE THIS IS IN THE 'wave' VARIABLE
    'zcoeffs', zCoeffs, 'measured pupil', pupilSize, ...
    'name', sprintf('human-%d', pupilSize),'spatial samples',nSamples); 
wvfP.calcpupilMM = pupilSize;
wvfP.refSizeOfFieldMM = refSizeOfFieldMM;

% Convert to siData format as well as wavefront object (feel free to change
% size of each sample)
[siPSFData, wvfP] = wvf2SiPsf(wvfP,'showBar',false,'nPSFSamples',nSamples,'umPerSample',1.1512); 

figure; 
subplot(1,2,1);
imagesc(wvfP.wavefrontaberrations{1}); 
axis square;
set(gca,'XTick',[]);
set(gca,'YTick',[]);
subplot(1,2,2);
imagesc(siPSFData.psf); 
axis square; 
colormap gray;
set(gca,'XTick',[]);
set(gca,'YTick',[]);