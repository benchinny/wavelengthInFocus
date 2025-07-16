%%

testa = imread('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/word_image_03.png');

 %%

testFFT = fftshift(fft2(squeeze(testa(:,:,1))));

xSize = floor(size(testa,2)/2);
ySize = floor(size(testa,1)/2);

[xx, yy] = meshgrid((1:size(testa,2))-xSize-1,(1:size(testa,1))-ySize-1);

dist = sqrt(xx.^2 + yy.^2);

A = [];
indCounter = [];
for j = 1:30
   ind = dist>(j-0.5) & dist<(j+0.5);
   A(j) = mean(abs(testFFT(ind)));
   indCounter(j) = sum(sum(ind));
end