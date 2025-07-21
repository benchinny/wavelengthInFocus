function wvInFocus = ARCwvInFocusConesMeanZbalance(subjNum,stimNum,wLMS,wQualityBalance)

wave = 380:4:780;
nFocus = length(wave);
foldernameCones = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImages/';

% USE THE SAME ORIGINAL (PRE-OPTICS) IMAGE EACH TIME--THIS ONE HAPPENS TO
% LIVE IN THE FOLDER FOR SUBJECT 10, BUT IT REALLY DOESN'T MATTER SINCE ALL
% SUBJECTS SAW THE SAME ON-SCREEN STIMULUS
fnameConeRspNoLCA = ['subj10block3stimulus1' 'focusInd1noLCA'];
absorptionsOrig = load([foldernameCones 'S10/' fnameConeRspNoLCA]);
absorptionsOrig = absorptionsOrig.absorptions;
coneImgOrig = sum(absorptionsOrig,3);
peakCorr = [];
LSbalance = [];

for i = 1:nFocus
    fnameConeRsp = ['subj' num2str(subjNum) 'stimulus' num2str(stimNum) 'focusInd' num2str(i)];
    load([foldernameCones 'S' num2str(subjNum) '/' fnameConeRsp]);
    % WEIGHT AND ADD CONE IMAGES
    coneImg = squeeze(absorptions(:,:,1).*wLMS(1) + absorptions(:,:,2).*wLMS(2) + absorptions(:,:,3).*wLMS(3));
    % PEAK CORRELATION OF 'QUALITY' IMAGE
    peakCorr(i) = max(max(normxcorr2(coneImgOrig,coneImg)));
    % CORRELATING L AND S IMAGES WITH ORIGINAL
    peakCorrL = max(max(normxcorr2(coneImgOrig,squeeze(absorptions(:,:,1)))));
    peakCorrS = max(max(normxcorr2(coneImgOrig,squeeze(absorptions(:,:,3)))));
    % CALCULATING DIFFERENCE IN IMAGE QUALITY BETWEEN L AND S IMAGES
    LSbalance(i) = peakCorrL-peakCorrS;
end

% WAVELENGTH TO FOCUS IF MAXIMIZING IMAGE QUALITY
[~,indPeakPeak] = max(peakCorr);
wvInFocusQuality = wave(indPeakPeak);
% WAVELENGTH TO FOCUS IF MAXIMIZING BALANCE
[~,indPeakBalance] = min(abs(LSbalance));
wvInFocusBalance = wave(indPeakBalance);

% CONVERT WAVELENGTHS TO FOCUS TO DIOPTERS, THEN WEIGHT AND ADD
wvInFocusDiopters = wQualityBalance(1)*humanWaveDefocusARC(550,wvInFocusQuality,99) + ...
                    wQualityBalance(2)*humanWaveDefocusARC(550,wvInFocusBalance,99);

% REVERSE TO OBTAIN WAVELENGTH THAT SHOULD BE FOCUSED
wvInFocus = humanWaveDefocusInvertARC(550,wvInFocusDiopters,99);

end