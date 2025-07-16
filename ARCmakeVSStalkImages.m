%%

im = imread('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/15532373947_3e6439e94d_o.jpg');

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/siPSFData.mat');

%%

scalePixFactor = 366.3;
psfInd = 51;

imBlurred = [];

imBlurred(:,:,1) = conv2(squeeze(im(:,:,1)),siPSFData.psf(:,:,psfInd));
imBlurred(:,:,2) = conv2(squeeze(im(:,:,2)),siPSFData.psf(:,:,psfInd));
imBlurred(:,:,3) = conv2(squeeze(im(:,:,3)),siPSFData.psf(:,:,psfInd));

imBlurred = imBlurred./scalePixFactor;

figure; 
imagesc(imBlurred);
set(gcf,'Position',[584 495 560*1 358*1]);

%%

clear all;

