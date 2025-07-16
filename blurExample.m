%%

x = linspace(0,1,60);
y = linspace(0,1,60);

theta = 0;
[xx, yy] = meshgrid(x,y);
xx2 = xx.*cos(theta)-yy.*sin(theta);
f = 1:20;

oi = oiCreateARC('human',580,0);
psf1 = fftshift(ifft2(oi.optics.OTF.OTF));
oi = oiCreateARC('human',540,0);
psf2 = fftshift(ifft2(oi.optics.OTF.OTF));
oi = oiCreateARC('human',510,0);
psf3 = fftshift(ifft2(oi.optics.OTF.OTF));

for i = 1:length(f)
    testImg = sin(2.*pi.*f(i).*xx2);
    testImgOpt1 = conv2(testImg,psf1);
    testImgOpt1 = testImgOpt1(31:90,31:90);
    c1(i) = sqrt(mean(testImgOpt1(:).^2));
    testImgOpt2 = conv2(testImg,psf2);
    testImgOpt2 = testImgOpt2(31:90,31:90);
    c2(i) = sqrt(mean(testImgOpt2(:).^2));    
    testImgOpt3 = conv2(testImg,psf3);
    testImgOpt3 = testImgOpt3(31:90,31:90);
    c3(i) = sqrt(mean(testImgOpt3(:).^2));        
end

figure; 
hold on; 
plot(log(f),c1,'LineWidth',1.5); 
plot(log(f),c2,'LineWidth',1.5);
plot(log(f),c3,'LineWidth',1.5);
axis square;
set(gca,'FontSize',15);
xlabel('Spatial Frequency (cyc/deg)');
ylabel('Contrast');
set(gca,'XTick',[]);
set(gca,'YTick',[]);

figure; 
hold on; 
plot(f,c1); 
plot(f,c2);
plot(f,c3);

