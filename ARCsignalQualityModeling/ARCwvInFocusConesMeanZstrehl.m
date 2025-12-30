function wvInFocus = ARCwvInFocusConesMeanZstrehl(subjNum,rgb,wLMS,dataPath)

% DEFINING WAVELENGTH
wave = 380:4:780;
nFocus = length(wave);

% PATH TO PSFs
foldernameCones = fullfile(dataPath,'data','psfs');
% PATH TO DISPLAY CALIBRATION FILES (TO GET PRIMARIES)
calPath = fullfile(dataPath,'data','helperFiles','BVAMS_calibration_files','Ben_calibration_July_6_2024');

% LOAD CALIBRATION FILES
load(fullfile(calPath,'redPrimaryJuly0624_initialPositionFocus3_100.mat'));
spd(:,1) = energy;
load(fullfile(calPath,'greenPrimaryJuly0624_initialPositionFocus3_100.mat'));
spd(:,2) = energy;
load(fullfile(calPath,'bluePrimaryJuly0624_initialPositionFocus3_100.mat'));
spd(:,3) = energy;

% GAMMA-CORRECT INPUT RGB VALUES
rgbGammaCorrected = rgb.^[2.5 2.7 2.3];
% CALCULATE STIMULUS ENERGY SPECTRUM FOR INPUT RGB VALUES
stimSPD = spd(:,1).*rgbGammaCorrected(1)+spd(:,2).*rgbGammaCorrected(2)+spd(:,3).*rgbGammaCorrected(3);

% LOAD CONE QUANTAL EFFICIENCIES
load(fullfile(dataPath,'data','helperFiles','sQE.mat'));
% LOAD TRANSMITTANCE VALUES
load(fullfile(dataPath,'data','helperFiles','transmittance.mat'));

% CALCULATE ENERGY THAT ACTUALLY GETS THROUGH THE LENS
stimSPD = reshape(stimSPD,[1 1 length(wave)]);
stimSPDtransmitted = bsxfun(@times,stimSPD,transmittance);

% INITIALIZE VECTOR FOR STORING NORMALIZED STREHL
strehlNorm = [];

for i = 1:nFocus % FOR EACH WAVELENGTH
    % LOAD POINT SPREAD FUNCTIONS
    fnameConeRsp = ['subj' num2str(subjNum) 'PSFfocusInd' num2str(i)];
    S = load(fullfile(foldernameCones,['S' num2str(subjNum)],fnameConeRsp));
    % SCALE EACH POINT SPREAD FUNCTION BY ENERGY GETTING PAST LENS
    psfIrradianceScaled = bsxfun(@times,S.psf,stimSPDtransmitted);
    % INITIALIZE MATRIX FOR STORING CONE IMAGES OF PSFS
    psfCone = [];
    for j = 1:size(sQE,2) % FOR EACH CONE TYPE
        % MULTIPLY EACH PSF BY THE CONE'S SENSITIVITY AT EACH WAVELENGTH
        sQEsingleCone = reshape(sQE(:,j),[1 1 length(wave)]);
        % THEN SUM TO GET THE 'CONE IMAGES' OF THE PSFS
        psfCone(:,:,j) = sum(bsxfun(@times,psfIrradianceScaled,sQEsingleCone),3);
    end
    % SUM THE CONE IMAGES. IF w_L = 0.72, w_M = 0.28, and w_S = 0, THAT
    % WILL BE EQUIVALENT TO THE LUMINANCE PSF
    psfMech = wLMS(1).*psfCone(:,:,1) + wLMS(2).*psfCone(:,:,2) + wLMS(3).*psfCone(:,:,3);
    % NORMALIZED STREHL IS JUST THE MAXIMUM OF THIS (WE DON'T DIVIDE BY THE
    % DIFFRACTION-LIMITED PSF PEAK BECAUSE IT MAKES NO DIFFERENCE FOR OUR
    % MODEL PREDICTION)
    strehlNorm(i) = max(psfMech(:));
end

% FIND THE WAVELENGTH YIELDING BEST NORMALIZED STREHL
[~,indPeakPeak] = max(strehlNorm);
wvInFocus = wave(indPeakPeak);

end