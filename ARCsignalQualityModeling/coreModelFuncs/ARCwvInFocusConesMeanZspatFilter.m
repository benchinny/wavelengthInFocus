function wvInFocus = ARCwvInFocusConesMeanZspatFilter(subjNum,stimNum,wLMS,dataPath)

% FOLDER WITH PRE-GENERATED CONE IMAGES
foldernameCones = fullfile(dataPath,'data','coneImages');
% FOLDER NAME WITH HELPER FILES
foldernameHelpers = fullfile(dataPath,'data','helperFiles');

wave = 380:4:780; % SAMPLED WAVELENGTHS IN MODEL
nFocus = length(wave); % NUMBER OF WAVELENGTHS MODELED

% LOAD THE ORIGINAL (PRE-OPTICS) IMAGE 
fnameConeRspNoLCA = 'coneImagesNoOptics';
absorptionsOrig = load(fullfile(foldernameHelpers,fnameConeRspNoLCA));
absorptionsOrig = absorptionsOrig.absorptions;
coneImgOrig = sum(absorptionsOrig,3);

% LOAD SPATIAL FILTER
load(fullfile(foldernameHelpers,'freqFilterARC.mat'));

% MASK FOR SIMULATING 'HOLE' IN S-CONE IMAGE. START OUT WITH SUPPORT.
[SconeMaskSupportXX, SconeMaskSupportYY] = meshgrid(-90:91,-90:91);
% START WITH A MASK WITH NO HOLE
SconeMask = ones(size(SconeMaskSupportXX));
% MAKE HOLE IN CENTER. THIS IS A HARD EDGE.
SconeMask(sqrt(SconeMaskSupportXX.^2 + SconeMaskSupportYY.^2)<22.5) = 0;
% NOW TO MAKE THE EDGE SOFT. WE ARE GOING TO CONVOLVE A BLUR KERNEL WITH
% THE ENTIRE MASK. START WITH THE SUPPORT OVER THE KERNEL.
[softEdgeSupportXX, softEdgeKernSupportYY] = meshgrid(linspace(-1,1,9));
% USE A 2D GAUSSIAN TO MAKE THE KERNEL
softEdgeKernCol = mvnpdf([softEdgeSupportXX(:) softEdgeKernSupportYY(:)],[0 0],[0.3^2 0; 0 0.3^2]); 
softEdgeKern = reshape(softEdgeKernCol,size(softEdgeSupportXX))./sum(softEdgeKernCol(:));
% CONVOLVE
SconeMaskSoft = conv2(SconeMask,softEdgeKern);
% REMOVE PADDING AT EDGES FROM THE RESULT OF THE CONVOLUTION
SconeMask = SconeMaskSoft(5:186,5:186);
% THE CONVOLUTION WILL CAUSE SOME VALUES TO BE BETWEEN 0 AND 1 AT THE
% EDGES. WE WANT THESE VALUES TO BE EXACTLY 1.
SconeMask(:,1:5) = 1;
SconeMask(:,178:182) = 1;
SconeMask(1:5,:) = 1;
SconeMask(178:182,:) = 1;

if wLMS(3) == 0 % IF NOT BLUE-YELLOW MECHANISM, NO MASK
    SconeMask = ones(size(SconeMask));
end
% OTHERWISE, APPLY S-CONE MASK
coneImgOrig = coneImgOrig.*SconeMask;

% FILTER CONE IMAGE ACCORDING TO FREQUENCIES KNOWN TO DRIVE ACCOMMODATION
coneImgOrigFFT = fftshift(fft2(coneImgOrig));
coneImgOrigFilteredFFT = coneImgOrigFFT.*freqFilterARC;
coneImgOrigFiltered = real(ifft2(ifftshift(coneImgOrigFilteredFFT)));

peakCorr = []; % INITIALIZE VECTOR FOR IMAGE QUALITY
for i = 1:nFocus % LOOP OVER WAVELENGTHS IN FOCUS
    fnameConeRsp = ['subj' num2str(subjNum) 'stimulus' num2str(stimNum) 'focusInd' num2str(i)];
    load(fullfile(foldernameCones,['S' num2str(subjNum)],fnameConeRsp));
    % APPLY CONE WEIGHTS AND SUM
    absorptions(:,:,1) = SconeMask.*absorptions(:,:,1).*wLMS(1);
    absorptions(:,:,2) = SconeMask.*absorptions(:,:,2).*wLMS(2);
    absorptions(:,:,3) = SconeMask.*absorptions(:,:,3).*wLMS(3);
    coneImg = sum(absorptions,3);
    % FILTER SIGNAL ACCORDING TO FREQUENCIES KNOWN TO DRIVE ACCOMMODATION
    coneImgFFT = fftshift(fft2(coneImg));
    coneImgFilteredFFT = coneImgFFT.*freqFilterARC;
    coneImgFiltered = real(ifft2(ifftshift(coneImgFilteredFFT)));
    % COMPUTE IMAGE QUALITY ACCORDING TO X-CORRELATION
    peakCorr(i) = max(max(abs(normxcorr2(coneImgFiltered,coneImgOrigFiltered))));
end

% IDENTIY WAVELENGTH THAT IS BEST TO PUT IN FOCUS
[~,indPeakPeak] = max(peakCorr);
wvInFocus = wave(indPeakPeak);

end