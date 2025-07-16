function wvInFocus = ARCwvInFocusConesMeanZstrehl(subjNum,rgb,wLMS)

% DEFINING WAVELENGTH
wave = 380:4:780;
nFocus = length(wave);
foldernameCones = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImages/';

% Setting up display properties
d = struct;

bUseBVAMScal = 1; % if using BVAMS calibration data

if strcmp(getenv('USER'),'benjaminchin')
    calPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/BVAMS_calibration_files/Ben_calibration_July_6_2024/';
     stimPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/';
     savePath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImages/S';
end

if strcmp(getenv('USER'),'benchin')
    calPath = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/BVAMS_calibration_files/Ben_calibration_July_6_2024/';
    stimPath = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/';
    savePath = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImages/S';
end

if bUseBVAMScal
    load([calPath 'redPrimaryJuly0624_initialPositionFocus3_100.mat']);
    d.spd(:,1) = energy;
    load([calPath 'greenPrimaryJuly0624_initialPositionFocus3_100.mat']);
    d.spd(:,2) = energy;
    load([calPath 'bluePrimaryJuly0624_initialPositionFocus3_100.mat']);
    d.spd(:,3) = energy;
end

% CALCULATE STIMULUS ENERGY SPECTRUM FOR INPUT RGB VALUES
rgbGammaCorrected = rgb.^[2.5 2.7 2.3];
stimSPD = d.spd(:,1).*rgbGammaCorrected(1)+d.spd(:,2).*rgbGammaCorrected(2)+d.spd(:,3).*rgbGammaCorrected(3);

% LOAD CONE QUANTAL EFFICIENCIES
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/helperFiles/sQE.mat');
% LOAD TRANSMITTANCE VALUES
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/helperFiles/transmittance.mat');

% CALCULATE ENERGY THAT ACTUALLY GETS THROUGH THE LENS
stimSPD = reshape(stimSPD,[1 1 length(wave)]);
stimSPDtransmitted = bsxfun(@times,stimSPD,transmittance);

% % Put the image center in (1, 1) and take the transform.
% imgFFT = fft2(fftshift(img));
% % Multiply the transformed otf and the image.
% % Then invert and put the image center in  the center of the matrix
% filteredIMG = abs(ifftshift(ifft2(otf .* imgFFT)));

strehlNorm = [];

for i = 1:nFocus % FOR EACH WAVELENGTH
    % LOAD POINT SPREAD FUNCTIONS
    fnameConeRsp = ['subj' num2str(subjNum) 'stimulus12focusInd' num2str(i) 'psf'];
    S = load([foldernameCones 'S' num2str(subjNum) '/' fnameConeRsp]);
    % SCALE EACH POINT SPREAD FUNCTION BY ENERGY GETTING PAST LENS
    psfIrradianceScaled = bsxfun(@times,S.psf,stimSPDtransmitted);
    psfCone = [];
    for j = 1:size(sQE,2) % FOR EACH CONE TYPE
        sQEsingleCone = reshape(sQE(:,j),[1 1 length(wave)]);
        psfCone(:,:,j) = sum(bsxfun(@times,psfIrradianceScaled,sQEsingleCone),3);
    end
    psfMech = wLMS(1).*psfCone(:,:,1) + wLMS(2).*psfCone(:,:,2) + wLMS(3).*psfCone(:,:,3);
    strehlNorm(i) = max(psfMech(:));
end

[~,indPeakPeak] = max(strehlNorm);
wvInFocus = wave(indPeakPeak);

end