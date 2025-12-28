function [wvInFocus, wave, peakCorr, rmsContrast] = ARCwvInFocusConesMeanZtestbed(subjNum,stimNum,wLMS,wv2examine)

% rgb00 = [0.3270         0    1.0000; ...
%          0.3270    0.3340    1.0000; ...
%          0.4320         0    1.0000; ...
%          0.4320    0.3340    1.0000; ...
%          0.4720         0    0.8150; ...
%          0.5690         0    0.5470; ...
%          0.5690         0    0.7400; ...
%          0.5690         0    1.0000; ...
%          0.5690    0.3340    0.5470; ...
%          0.5690    0.3340    0.7400; ...
%          0.5690    0.3340    1.0000; ...
%          0.5690    0.4320    1.0000];

wave = 380:4:780;
nFocus = length(wave);
foldernameCones = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImages/';

% USE THE SAME ORIGINAL (PRE-OPTICS) IMAGE EACH TIME--THIS ONE HAPPENS TO
% LIVE IN THE FOLDER FOR SUBJECT 10, BUT IT REALLY DOESN'T MATTER SINCE ALL
% SUBJECTS SAW THE SAME ON-SCREEN STIMULUS
fnameConeRspNoLCA = ['subj10block3stimulus1' 'focusInd1noLCA'];
absorptionsOrig = load([foldernameCones 'S10/' fnameConeRspNoLCA]);
absorptionsOrig = absorptionsOrig.absorptions;
% JUST ADD THEM ALL UP--THERE IS NO DEFOCUS ANYWAY
coneImgOrig = sum(absorptionsOrig,3);
peakCorr = [];
rmsContrast = [];

for i = 1:nFocus
    fnameConeRsp = ['subj' num2str(subjNum) 'stimulus' num2str(stimNum) 'focusInd' num2str(i)];
    load([foldernameCones 'S' num2str(subjNum) '/' fnameConeRsp]);
    absorptions(:,:,1) = absorptions(:,:,1).*wLMS(1);
    absorptions(:,:,2) = absorptions(:,:,2).*wLMS(2);
    absorptions(:,:,3) = absorptions(:,:,3).*wLMS(3);
    coneImg = sum(absorptions,3);
    peakCorr(i) = max(max(normxcorr2(coneImgOrig,coneImg)));
    % TESTING OUT MEASURES OF LUMINANCE / CONTRAST
    coneImgThresh = coneImg;
    % rmsContrast(i) = sqrt(mean((coneImgThresh(:)-mean(coneImgThresh(:))).^2));
    % rmsContrast(i) = max(coneImgThresh(:));
    rmsContrast(i) = sum(coneImgThresh(:));
    % PLOT CONE IMAGES FOR SPECIFIED WAVELENGTHS IN FOCUS
    if ismember(wave(i),wv2examine)
        figure;
        set(gcf,'Position',[377 343 701 603]);
        subplot(2,2,1);
        imagesc(squeeze(absorptions(:,:,1)));
        colormap gray;
        colorbar;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        title('L-cones');
        subplot(2,2,2);
        imagesc(squeeze(absorptions(:,:,2)));
        colormap gray;
        colorbar;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]); 
        title('M-cones');
        subplot(2,2,3);
        imagesc(squeeze(absorptions(:,:,3)));
        colormap gray;
        colorbar;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]); 
        title('S-cones');
        subplot(2,2,4);
        imagesc(coneImg);
        colormap gray;
        colorbar;
        set(gca,'XTick',[]);
        set(gca,'YTick',[]); 
        title(['Wavelength in focus = ' num2str(wave(i))]);
    end
end

% READ OFF THE PEAK
[~,indPeakPeak] = max(peakCorr);
wvInFocus = wave(indPeakPeak);

end