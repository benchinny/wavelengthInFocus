function wvInFocus = ARCwvInFocusConesMeanZdeltaPass(subjNum,rgb,wLMS,dataPath)

% SPATIAL FREQUENCY TO CHECK THE MTF AT
frqCpd = 2;

if mod(frqCpd,1)~=0
    error('ARCwvInFocusConesMeanZdeltaPass: frqCpd has to be integer for now');
end

% DEFINING WAVELENGTH
wave = 380:4:780;
nFocus = length(wave);
% PATH TO PSFs
foldernameCones = fullfile(dataPath,'data','psfs');
% PATH TO DISPLAY CALIBRATION FILES (TO GET PRIMARIES)
calPath = fullfile(dataPath,'data','helperFiles','BVAMS_calibration_files','Ben_calibration_July_6_2024');

% Setting up display properties
d = struct;
% LOADING DISPLAY CALIBRATION DATA
load(fullfile(calPath,'redPrimaryJuly0624_initialPositionFocus3_100.mat'));
d.spd(:,1) = energy;
load(fullfile(calPath,'greenPrimaryJuly0624_initialPositionFocus3_100.mat'));
d.spd(:,2) = energy;
load(fullfile(calPath,'bluePrimaryJuly0624_initialPositionFocus3_100.mat'));
d.spd(:,3) = energy;

% CALCULATE STIMULUS ENERGY SPECTRUM FOR INPUT RGB VALUES
rgbGammaCorrected = rgb.^[2.5 2.7 2.3];
stimSPD = d.spd(:,1).*rgbGammaCorrected(1)+d.spd(:,2).*rgbGammaCorrected(2)+d.spd(:,3).*rgbGammaCorrected(3);

% LOAD CONE QUANTAL EFFICIENCIES
load(fullfile(dataPath,'data','helperFiles','sQE.mat'));
% LOAD TRANSMITTANCE VALUES
load(fullfile(dataPath,'data','helperFiles','transmittance.mat'));

% CALCULATE ENERGY THAT ACTUALLY GETS THROUGH THE LENS
stimSPD = reshape(stimSPD,[1 1 length(wave)]);
stimSPDtransmitted = bsxfun(@times,stimSPD,transmittance);

% INITIALIZE VECTOR FOR STORING CONTRAST AT SPATIAL FREQ OF INTEREST
contrastAtFrqCpd = [];
% NEED TO SCALE BY SIZE OF ORIGINAL PSF GRID--'FREQUENCY IN CYLES PER
% STIMULUS SPACE'
frqCpStimSpace = frqCpd.*1.219;
% MESHGRID OVER FREQUENCIES OF OTF
fx = -29:30;
fy = -29:30;
[fxx, fyy] = meshgrid(fx,fy);
% GO IN A CIRCLE IN FREQUENCY SPACE AT THE FREQUENCY OF INTEREST
angle2interp = 0:360;
% GET THE CORRESPONDING X AND Y COORDINATES
xcoord = frqCpStimSpace.*cosd(angle2interp);
ycoord = frqCpStimSpace.*sind(angle2interp);

for i = 1:nFocus % FOR EACH WAVELENGTH
    % LOAD POINT SPREAD FUNCTIONS
    fnameConeRsp = ['subj' num2str(subjNum) 'PSFfocusInd' num2str(i)];
    S = load(fullfile(foldernameCones,['S' num2str(subjNum)],fnameConeRsp));
    % SCALE EACH POINT SPREAD FUNCTION BY ENERGY GETTING PAST LENS
    mtfIrradianceScaled = bsxfun(@times,abs(S.otf),stimSPDtransmitted);
    mtfCone = [];
    for j = 1:size(sQE,2) % FOR EACH CONE TYPE
        % MULTIPLY BY CONE FUNDAMENTALS
        sQEsingleCone = reshape(sQE(:,j),[1 1 length(wave)]);
        mtfCone(:,:,j) = sum(bsxfun(@times,mtfIrradianceScaled,sQEsingleCone),3);
    end
    % GET THE 'MECHANISM PSF' (E.G. THE LUMINANCE PSF)
    mtfMech = wLMS(1).*mtfCone(:,:,1) + wLMS(2).*mtfCone(:,:,2) + wLMS(3).*mtfCone(:,:,3);

    % MTF IS TYPICALLY NORMALIZED
    energy0(i) = mtfMech(size(mtfMech,1)/2,size(mtfMech,2)/2);
    % COMPUTE MEAN CONTRAST AT DISTANCE FRQCPD FROM ORIGIN
    contrastAtFrqCpd(i) = mean(interp2(fxx,fyy,mtfMech,xcoord,ycoord))./mtfMech(size(mtfMech,1)/2,size(mtfMech,2)/2);
end

% REMOVE VALUES FOR WHICH AMPLITUDE AT 0 IS TOO SMALL--CAUSES INSTABILITIES
energy0threshold = 0.0004;
contrastAtFrqCpdInRange = contrastAtFrqCpd(energy0>energy0threshold);
waveInRange = wave(energy0>energy0threshold);
% FIND WAVELENGTH THAT YIELDS PEAK CONTRAST IN RANGE OF INTEREST
[~,indPeakPeak] = max(contrastAtFrqCpdInRange);
wvInFocus = waveInRange(indPeakPeak);

end