%%

% im = imread('/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/My Drive/K99/Figures/MalteseCross2.png');
im = AFCwordStim('car',[500 500],[70 70; 160 70; 260 70],'blue',200);

rScale = 1;
bScale = 1.00;

imR = double(im);
imR(:,:,2) = 0;
imR(:,:,3) = 0;
imR = round(imR.*rScale);

imB = double(im);
imB(:,:,1) = 0;
imB(:,:,2) = 0;
imB = round(imB.*bScale);
figure; imagesc((imB./255).^(1/2.2)); set(gca,'XTick',[]); set(gca,'YTick',[]); clim([0 1]);
%%
[xx, yy] = meshgrid(-1:0.025:1,-1:0.025:1);
blurKernel = mvnpdf([xx(:) yy(:)],[0 0],[0.1 0; 0 0.1]);
blurKernel = blurKernel./sum(blurKernel(:));
blurKernel = reshape(blurKernel,size(xx));
imBlurB = conv2(squeeze(imB(:,:,3)),blurKernel);
imBlurB = imBlurB(41:540,41:540);
imBlur(:,:,1) = zeros(size(imBlurB));
imBlur(:,:,2) = zeros(size(imBlurB));
imBlur(:,:,3) = imBlurB;

% imBlur2 = zeros(size(imBlur));
% imBlur2(:,:,1) = imBlurB;

blurKernel2 = mvnpdf([xx(:) yy(:)],[0 0],[0.1 0; 0 0.1]);
blurKernel2 = blurKernel2./sum(blurKernel2(:));
blurKernel2 = reshape(blurKernel2,size(xx));
imBlurB2 = conv2(squeeze(imR(:,:,1)),blurKernel2);
imBlurB2 = imBlurB2(41:540,41:540);
imBlur2(:,:,3) = zeros(size(imBlurB2));
imBlur2(:,:,2) = zeros(size(imBlurB2));
imBlur2(:,:,1) = imBlurB2;
% 
imG = zeros(size(imBlur));
imG(:,:,2) = im(:,:,2);

figure; imagesc((imBlur./255).^(1/2.2)); set(gca,'XTick',[]); set(gca,'YTick',[]); clim([0 1]);

figure; imagesc(((1*imBlur+imBlur2.*1)./255).^(1/2.2)); set(gca,'XTick',[]); set(gca,'YTick',[]); clim([0 1]); set(gcf,'Position',[680 503 560 490])

%%

x = 0:0.01:6.5;
y = imresize([3 4],[1 length(x)],'nearest');
figure;
set(gcf,'Position',[387 486 1103 420]);
subplot(1,2,1);
hold on;
plot(x(1:round(length(x)/2)),y(1:round(length(x)/2)),'-','LineWidth',2,'Color',[0.62 0 0]);
plot(x((round(length(x)/2)+1):end),y((round(length(x)/2)+1):end),'-','LineWidth',2,'Color',[0.62 0 1]);
xlim([0 6.5]);
ylim([0 5]);
set(gca,'FontSize',20);
xlabel('Time(s)'); 
ylabel('Accommodative demand (D)');
y = imresize([3 2],[1 length(x)],'nearest');
subplot(1,2,2);
hold on;
plot(x(1:round(length(x)/2)),y(1:round(length(x)/2)),'-','LineWidth',2,'Color',[0.62 0 0]);
plot(x((round(length(x)/2)+1):end),y((round(length(x)/2)+1):end),'-','LineWidth',2,'Color',[0.62 0 1]);
xlim([0 6.5]);
ylim([0 5]);
set(gca,'FontSize',20);

