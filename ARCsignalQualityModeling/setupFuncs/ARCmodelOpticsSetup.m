function [oi, siPSFData]= ARCmodelOpticsSetup(subjNum,zCoeffs,calcWavelengths,measuredWavelength,PupilSize,spatialSamplesXY,defocusSet)

% This function generates an optics struct for ISETBIO to model the retinal
% image in the accommodation and acuity tasks (Exp 1 and Exp 2
% respectively). WARNING: the settings in this function are specific to the
% BVAMS apparatus used at UC Berkeley.
%
% subjNum: subject number (valid numbers: 1, 3, 5, 10, 16, 17, 18, 20)
% zCoeffs: Zernike coefficients
% calcWavelengths: what wavelengths to model the retinal image at
% measured wavelength: the wavelength in focus
% PupilSize: pupil size
% spatialSamplesXY: size of the stimulus in X and Y (not rows and columns!)
% defocusSet: value to manually set the defocus coefficient at

% MICROMETERS PER SAMPLE. FOR THE CHIN ET AL. (2026) MODEL, THIS WILL YIELD
% 0.23 ARCMINS PER PIXEL. NOTE THAT THIS VALUE IS SPECIFIC TO THE CODE FOR
% CHIN ET AL.
umPerSample = 1.1512;

% CREATE ISETBIO WAVEFRONT OBJECT
wvfP = wvfCreate('calc wavelengths', calcWavelengths, ...
    'measured wavelength', measuredWavelength, ...
    'zcoeffs', zCoeffs, 'measured pupil', PupilSize, ...
    'name', sprintf('human-%d', PupilSize),'spatial samples',spatialSamplesXY(1));
% MAKE SURE THE calcpupilMM FIELD MATCHES THE ACTUAL PUPIL SIZE IN
% THE EXPERIMENT
wvfP.calcpupilMM = PupilSize;
% SET CUSTOM LCA FUNCTION PER SUBJECT--ISETBIO WANTS IT TO BE SET
% IN A PARTICULAR FORMAT
if subjNum==1
    % CALCULATE MINIMUM AND MAXIMUM DEFOCUS FROM LCA FOR A
    % PARTICULAR SUBJECT FOR THE RANGE OF WAVELENGTHS TO ANALYZE.
    % THIS WILL ENSURE THAT THE MESH OVER THE PUPIL FUNCTION SPANS
    % A SUFFICIENT RANGE (refSizeOfFieldMM). IF DEFOCUS IS LARGE
    % ENOUGH, THE RANGE NEEDS TO BE REDUCED.
    defocusFromLCA = max(abs([humanWaveDefocusS1(measuredWavelength,min(calcWavelengths)) ...
                              humanWaveDefocusS1(measuredWavelength,max(calcWavelengths))]+defocusSet));
    wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS1);
elseif subjNum==3
    defocusFromLCA = max(abs([humanWaveDefocusS3(measuredWavelength,min(calcWavelengths)) ...
                              humanWaveDefocusS3(measuredWavelength,max(calcWavelengths))]+defocusSet));  
    wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS3);
elseif subjNum==5
    defocusFromLCA = max(abs([humanWaveDefocusS5(measuredWavelength,min(calcWavelengths)) ...
                              humanWaveDefocusS5(measuredWavelength,max(calcWavelengths))]+defocusSet));  
    wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS5); 
elseif subjNum==10
    defocusFromLCA = max(abs([humanWaveDefocusS10(measuredWavelength,min(calcWavelengths)) ...
                              humanWaveDefocusS10(measuredWavelength,max(calcWavelengths))]+defocusSet));  
    wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS10); 
elseif subjNum==16
    defocusFromLCA = max(abs([humanWaveDefocusS16(measuredWavelength,min(calcWavelengths)) ...
                              humanWaveDefocusS16(measuredWavelength,max(calcWavelengths))]+defocusSet));  
    wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS16); 
elseif subjNum==17
    defocusFromLCA = max(abs([humanWaveDefocusS17(measuredWavelength,min(calcWavelengths)) ...
                              humanWaveDefocusS17(measuredWavelength,max(calcWavelengths))]+defocusSet));  
    wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS17); 
elseif subjNum==18
    defocusFromLCA = max(abs([humanWaveDefocusS18(measuredWavelength,min(calcWavelengths)) ...
                              humanWaveDefocusS18(measuredWavelength,max(calcWavelengths))]+defocusSet));  
    wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS18); 
elseif subjNum==20
    defocusFromLCA = max(abs([humanWaveDefocusS20(measuredWavelength,min(calcWavelengths)) ...
                              humanWaveDefocusS20(measuredWavelength,max(calcWavelengths))]+defocusSet));  
    wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS20); 
else
    error('Subject number with no LCA function?');
end

% IF DEFOCUS IS LARGE ENOUGH, THE AREA OF THE PUPIL THE PSF IS
% CALCULATED FROM NEEDS TO BE REDUCED, OR YOU WILL END UP WITH A
% DEGENERATE PSF
if defocusFromLCA<1
    wvfP.refSizeOfFieldMM = 12;
else
    wvfP.refSizeOfFieldMM = 6;
end
% SET THE COEFFICIENT ON DEFOCUS EXPLICITLY
wvfP = wvfSet(wvfP, 'zcoeff', defocusSet, 'defocus');

% MAKE POINT-SPREAD FUNCTION (siPSFData) AND WAVEFRONT STRUCT
[siPSFData, wvfP] = wvf2SiPsfARC(wvfP,'showBar',false,'nPSFSamples',spatialSamplesXY(1),'umPerSample',umPerSample); 
oi = wvf2oi(wvfP); % CONVERT WAVEFRONT STRUCT TO OPTICS OBJECT
% NEED TO REMOVE PADDED ZEROS FROM PSF TO MAKE SAME SIZE AS
% STIMULUS IMAGE. THE LINES BELOW IDENTIFY THE 'GOOD INDICES', I.E.
% THE INDICES THAT AREN'T THE PADDED ZEROS TO BE REMOVED
paddingXCpsf = round((size(siPSFData.psf,2)-spatialSamplesXY(1))/2);
paddingYRpsf = round((size(siPSFData.psf,1)-spatialSamplesXY(2))/2);
indNotPadded = {(paddingYRpsf+1):(size(siPSFData.psf,1)-paddingYRpsf) ...
                (paddingXCpsf+1):(size(siPSFData.psf,2)-paddingXCpsf)};
% REPLACE OTF FIELD WITH A NEW OTF CALCULATED FROM THE PSFS WE
% JUST CALCULATED. COULDN'T FIGURE OUT HOW TO MAKE ISETBIO DO
% THIS AUTOMATICALLY, SO I DID IT MANUALLY (AND CHECKED THE
% OUTPUT!)
oi.optics.OTF = []; % INITIALIZE ARRAY FOR STORING OTFS
for j = 1:size(siPSFData.psf,3) % LOOP OVER WAVELENGTHS
    % NOTE THAT WE APPLY fftshift TO THE PSF SO THAT ITS CENTER IS
    % IN INDEX (1,1) OF THE IMAGE (TOP LEFT). fft2 EXPECTS THE
    % SIGNAL ORIGIN TO BE IN THIS LOCATION.
    oi.optics.OTF.OTF(:,:,j) = fft2(fftshift(squeeze(siPSFData.psf(indNotPadded{1},indNotPadded{2},j))));
end

end