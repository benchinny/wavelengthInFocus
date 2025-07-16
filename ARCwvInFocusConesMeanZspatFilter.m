function wvInFocus = ARCwvInFocusConesMeanZspatFilter(subjNum,stimNum,wLMS)

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

% LOAD SPATIAL FILTER
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/modelParams/freqFilterARC.mat');

coneImgOrigFFT = fftshift(fft2(coneImgOrig));
coneImgOrigFilteredFFT = coneImgOrigFFT.*freqFilterARC;
coneImgOrigFiltered = real(ifft2(ifftshift(coneImgOrigFilteredFFT)));

coneImgOrigFFT2 = fft2(fftshift(coneImgOrig));
freqFilterARCfft2 = fftshift(freqFilterARC);
coneImgOrigFilteredFFT2 = coneImgOrigFFT2.*freqFilterARCfft2;
coneImgOrigFiltered2 = real(ifftshift(ifft2(coneImgOrigFilteredFFT2)));

% % Put the image center in (1, 1) and take the transform.
% imgFFT = fft2(fftshift(img));
% % Multiply the transformed otf and the image.
% % Then invert and put the image center in  the center of the matrix
% filteredIMG = abs(ifftshift(ifft2(otf .* imgFFT)));

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

    % peakCorr(i) = max(max(abs(normxcorr2(coneImgFiltered,coneImgOrigFiltered))));
    peakCorr(i) = max(max(normxcorr2(coneImgFiltered,coneImgOrigFiltered)));
end

[~,indPeakPeak] = max(peakCorr);
wvInFocus = wave(indPeakPeak);

end