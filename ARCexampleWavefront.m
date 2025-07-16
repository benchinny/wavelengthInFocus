%%

PARAMS.PixelDimension = 512;% size of pupil aperture field in pixels (this defines the resolution of the calculation)
PARAMS.PupilSize = 7; %default values - will be replaced depending on choices below
PARAMS.PupilFieldSize =42; %default values - will be replaced depending on choices below
PARAMS.PupilFitSize = 7; %default values - will be replaced depending on choices below
PARAMS.ImagingWavelength = 0.55;% imaging wavelength in microns
PARAMS.WavefrontResolution = 53;% increase to enhance the display of the wavefront (doesn't affect calculation)

saveWF = 0;  %set to zero to skip this part of the code
savePSF = 0; %set to zero to skip this part of the code
runCONV = 1; %set to zero to skip this part of the code
runMTFPTF = 1; %set to zero to skip this part of the code
saveMTF = 1;
ProcessWholeMovie =0; %specific to reading in FIAT zernike lists (see below)

[fname pname] = uigetfile('*.CSV;*.csv', 'Select FIAT video-processing output CSV file');% open CSV analysis file from the ForsightV6 system
fpresize = size(fname,2);
fext = fname(fpresize-2:end);
fnametemp = fname(1:fpresize-4);

TimeStamp = datetime;
TimeStamp.Format = 'yy-MM-dd-hh-mm';
TimeStamp = char(TimeStamp);

ZernikeTable = readtable([pname '/' fname]); % read in the CSV file from the ForsightV6 system

FrameStart = 1; %first frame for analysis
FrameEnd = 1; %last frame for analysis
NumFrames = FrameEnd - FrameStart + 1;

NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 

c=zeros(1,65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 

%set these parameter based on what is contained in the *.zer file
PARAMS.PupilSize=table2array(ZernikeTable(FrameStart,5)); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
PARAMS.PupilFitSize=table2array(ZernikeTable(FrameStart,5)); 
PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
fprintf('\nThe current *.zer filename is %s\n',fname);
PARAMS.FileName = fname;
c(3:NumCoeffs)=table2array(ZernikeTable(FrameStart,11:width(ZernikeTable)));

%% Create wavefront structure with reasonable parameters.
pupilMM = 4;
% zCoeffs = wvfLoadThibosVirtualEyes(pupilMM);
% zCoeffs = [0 c(1:end-1)];
zCoeffs = zeros([1 65]);
zCoeffs(1) = 0;
zCoeffs(2) = 0;
zCoeffs(3) = 0;
zCoeffs(4) = 0.1;
zCoeffs(5) = 0.2887*0;
wave = [780];
wvfP = wvfCreate('calc wavelengths', wave, ...xl
    'measured wavelength', 780, ... % THIS REALLY MEANS 'REFERENCE' WAVELENGTH
    'zcoeffs', zCoeffs, 'measured pupil', pupilMM, ...
    'name', sprintf('human-%d', pupilMM),'spatial samples',320);
wvfP.calcpupilMM = pupilMM;
wvfP.refSizeOfFieldMM = 11;

% Set a little defocus, just to make the PSF a bit more interesting
% wvfP = wvfSet(wvfP, 'zcoeff', 0, 'defocus');

% Convert to siData format and save.  201 is the number of default 
% samples in the wvfP object, and we need to match that here.
[siPSFData, wvfP] = wvf2SiPsf(wvfP,'showBar',false,'nPSFSamples',320,'umPerSample',1.1512);

figure; 
set(gcf,'Position',[318 486 928 420]);
for i = 1:length(wave)
   subplot(1,2,1);
   imagesc(flipud(siPSFData.psf(:,:,i))); 
   colormap gray; 
   axis square
   title(['wave = ' num2str(wave(i)) 'mm']);
   subplot(1,2,2);
   plot(siPSFData.psf(161,:,1)); 
   colormap gray; 
   axis square
   title(['wave = ' num2str(wave(i)) 'mm']);   
   pause;
end

