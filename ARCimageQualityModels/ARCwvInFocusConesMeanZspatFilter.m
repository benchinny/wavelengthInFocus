function wvInFocus = ARCwvInFocusConesMeanZspatFilter(subjNum,stimNum,wLMS,dataPath)

if ispc
    slash = '\';
else
    slash = '/';
end
foldernameCones = [dataPath 'data' slash 'coneImages' slash];

wave = 380:4:780;
nFocus = length(wave);

% USE THE SAME ORIGINAL (PRE-OPTICS) IMAGE EACH TIME--THIS ONE HAPPENS TO
% LIVE IN THE FOLDER FOR SUBJECT 10, BUT IT REALLY DOESN'T MATTER SINCE ALL
% SUBJECTS SAW THE SAME ON-SCREEN STIMULUS
fnameConeRspNoLCA = ['subj10block3stimulus1' 'focusInd1noLCA'];
absorptionsOrig = load([foldernameCones 'S10' slash fnameConeRspNoLCA]);
absorptionsOrig = absorptionsOrig.absorptions;
coneImgOrig = sum(absorptionsOrig,3);

% LOAD SPATIAL FILTER
load([dataPath 'data' slash 'modelParams' slash 'freqFilterARC.mat']);

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

[SconeMaskSupportXX, SconeMaskSupportYY] = meshgrid(-90:91,-90:91);
SconeMask = ones(size(SconeMaskSupportXX));
SconeMask(sqrt(SconeMaskSupportXX.^2 + SconeMaskSupportYY.^2)<22.5) = 0;
[softEdgeSupportXX, softEdgeKernSupportYY] = meshgrid(linspace(-1,1,9));
softEdgeKernCol = mvnpdf([softEdgeSupportXX(:) softEdgeKernSupportYY(:)],[0 0],[0.3^2 0; 0 0.3^2]); 
softEdgeKern = reshape(softEdgeKernCol,size(softEdgeSupportXX))./sum(softEdgeKernCol(:));
SconeMaskSoft = conv2(SconeMask,softEdgeKern);
SconeMask = SconeMaskSoft(5:186,5:186);
SconeMask(:,1:5) = 1;
SconeMask(:,178:182) = 1;
SconeMask(1:5,:) = 1;
SconeMask(178:182,:) = 1;

peakCorr = [];
for i = 1:nFocus
    fnameConeRsp = ['subj' num2str(subjNum) 'stimulus' num2str(stimNum) 'focusInd' num2str(i)];
    load([foldernameCones 'S' num2str(subjNum) slash fnameConeRsp]);
    absorptions(:,:,1) = absorptions(:,:,1).*wLMS(1);
    absorptions(:,:,2) = absorptions(:,:,2).*wLMS(2);
    absorptions(:,:,3) = SconeMask.*absorptions(:,:,3).*wLMS(3);
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