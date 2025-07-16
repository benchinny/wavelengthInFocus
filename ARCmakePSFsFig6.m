%%

egPSF1tmp = fftshift(ifft2(squeeze(oi.optics.OTF.OTF(:,:,61))));
egPSF1 = [];
egPSF1(:,:,1) = egPSF1tmp.*255/2;
egPSF1(:,:,2) = egPSF1tmp.*0;
egPSF1(:,:,3) = egPSF1tmp.*0;
figure; imagesc(egPSF1); axis square; 
set(gca,'XTick',[]);
set(gca,'YTick',[]);

egPSF2tmp = fftshift(ifft2(squeeze(oi.optics.OTF.OTF(:,:,51))));
egPSF2 = [];
egPSF2(:,:,1) = egPSF2tmp.*255/2;
egPSF2(:,:,2) = egPSF2tmp.*128/2;
egPSF2(:,:,3) = egPSF2tmp.*0;
figure; imagesc(egPSF2); axis square; 
set(gca,'XTick',[]);
set(gca,'YTick',[]);

egPSF3tmp = fftshift(ifft2(squeeze(oi.optics.OTF.OTF(:,:,41))));
egPSF3 = [];
egPSF3(:,:,1) = egPSF3tmp.*0;
egPSF3(:,:,2) = egPSF3tmp.*255/2;
egPSF3(:,:,3) = egPSF3tmp.*0;
figure; imagesc(egPSF3); axis square; 
set(gca,'XTick',[]);
set(gca,'YTick',[]);

egPSF4tmp = fftshift(ifft2(squeeze(oi.optics.OTF.OTF(:,:,31))));
egPSF4 = [];
egPSF4(:,:,1) = egPSF4tmp.*0;
egPSF4(:,:,2) = egPSF4tmp.*80;
egPSF4(:,:,3) = egPSF4tmp.*150;
figure; imagesc(egPSF4); axis square; 
set(gca,'XTick',[]);
set(gca,'YTick',[]);

egPSF5tmp = fftshift(ifft2(squeeze(oi.optics.OTF.OTF(:,:,21))));
egPSF5 = [];
egPSF5(:,:,1) = egPSF5tmp.*0;
egPSF5(:,:,2) = egPSF5tmp.*100;
egPSF5(:,:,3) = egPSF5tmp.*510;
figure; imagesc(egPSF5); axis square; 
set(gca,'XTick',[]);
set(gca,'YTick',[]);

%%

egPSF3tmp = fftshift(ifft2(squeeze(oi.optics.OTF.OTF(:,:,41))));
egPSF4tmp = fftshift(ifft2(squeeze(oi.optics.OTF.OTF(:,:,31))));
egPSF5tmp = fftshift(ifft2(squeeze(oi.optics.OTF.OTF(:,:,21))));

