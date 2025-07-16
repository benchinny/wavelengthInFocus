function wvInFocus = ARCwvInFocusConesMeanZnoAbbSpatFilter(subjNum,stimNum,wLMS)

% function wvInFocus = ARCwvInFocusConesMeanZnoAbbSpatFilter(subjNum,stimNum,wLMS)
%
% determining wavelength in focus for spatially filtered stimulus, with no
% aberrations and standard LCA

% INPUT VARIABLE stimNum INDICES FOR TESTING: 
% 1 : 0.3270         0    1.0000
% 2 : 0.3270    0.3340    1.0000
% 3 : 0.4320         0    1.0000
% 4 : 0.4320    0.3340    1.0000
% 5 : 0.4720         0    0.8150
% 6 : 0.5690         0    0.5470
% 7 : 0.5690         0    0.7400
% 8 : 0.5690         0    1.0000
% 9 : 0.5690    0.3340    0.5470
% 10: 0.5690    0.3340    0.7400
% 11: 0.5690    0.3340    1.0000
% 12: 0.5690    0.4320    1.0000

wave = 380:4:780;
nFocus = length(wave);
foldernameCones = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImagesNoAbb/';

% USE THE SAME ORIGINAL (PRE-OPTICS) IMAGE EACH TIME--THIS ONE HAPPENS TO
% LIVE IN THE FOLDER FOR SUBJECT 10, BUT IT REALLY DOESN'T MATTER SINCE ALL
% SUBJECTS SAW THE SAME ON-SCREEN STIMULUS
fnameConeRspNoLCA = ['subj10block3stimulus1' 'focusInd1noLCA'];
absorptionsOrig = load([foldernameCones 'S10/' fnameConeRspNoLCA]);
absorptionsOrig = absorptionsOrig.absorptions;
coneImgOrig = sum(absorptionsOrig,3);

% LOAD SPATIAL FILTER
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/modelParams/freqFilterARC.mat');

coneImgOrigFFT = fftshift(fft2(coneImgOrig));
coneImgOrigFilteredFFT = coneImgOrigFFT.*freqFilterARC;
coneImgOrigFiltered = real(ifft2(ifftshift(coneImgOrigFilteredFFT)));

peakCorr = [];

for i = 1:nFocus
    fnameConeRsp = ['subj' num2str(subjNum) 'stimulus' num2str(stimNum) 'focusInd' num2str(i)];
    load([foldernameCones 'S' num2str(subjNum) '/' fnameConeRsp]);
    absorptions(:,:,1) = absorptions(:,:,1).*wLMS(1);
    absorptions(:,:,2) = absorptions(:,:,2).*wLMS(2);
    absorptions(:,:,3) = absorptions(:,:,3).*wLMS(3);
    coneImg = sum(absorptions,3);

    coneImgFFT = fftshift(fft2(coneImg));
    coneImgFilteredFFT = coneImgFFT.*freqFilterARC;
    coneImgFiltered = real(ifft2(ifftshift(coneImgFilteredFFT)));
    peakCorr(i) = max(max(normxcorr2(coneImgOrigFiltered,coneImgFiltered)));
end

[~,indPeakPeak] = max(peakCorr);
wvInFocus = wave(indPeakPeak);

end