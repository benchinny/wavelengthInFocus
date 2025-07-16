%%

figure; 
set(gcf,'Position',[289 428 1056 420]);
subplot(1,3,1);
plot(d.wave,d.spd(:,1)'.*T_sensorXYZ(2,:).*0.2347.*0.25,'r','LineWidth',1.5); hold on;
plot(d.wave,d.spd(:,2)'.*T_sensorXYZ(2,:).*0.1037.*0.0,'g','LineWidth',1.5);
plot(d.wave,d.spd(:,3)'.*T_sensorXYZ(2,:),'b','LineWidth',1.5);
plot(d.wave,d.spd(:,1)'.*T_sensorXYZ(2,:).*0.2347.*0.25 + ...
     d.spd(:,2)'.*T_sensorXYZ(2,:).*0.1037.*0.0 + ...
     d.spd(:,3)'.*T_sensorXYZ(2,:),'-k','LineWidth',1.5);
axis square;
formatFigure('Wavelength (\lambda)','Luminance');
% subplot(1,3,2);
% imagesc(I);
% set(gca,'XTick',[]);
% set(gca,'YTick',[]);
% axis square;
% set(gca,'FontSize',15);
% title('Original');
subplot(1,3,3);
plot(s.spectrum.wave,squeeze(s.data.photons(142,160,:)),'-k','LineWidth',1);
formatFigure('Wavelength (\lambda)','Photons');
axis square;

%%

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/BVAMS_calibration_files/Ben_calibration_July_6_2024/redPrimaryJuly0624_initialPositionFocus3_100.mat');
spectralAxisInterp = 380:780;
energyInterp = interp1(spectralAxis,energy,spectralAxisInterp,'spline');
figure; plot(spectralAxisInterp,energyInterp,'r-');
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/BVAMS_calibration_files/Ben_calibration_July_6_2024/greenPrimaryJuly0624_initialPositionFocus3_100.mat');
spectralAxisInterp = 380:780;
energyInterp = interp1(spectralAxis,energy,spectralAxisInterp,'spline');
figure; plot(spectralAxisInterp,energyInterp,'g-');
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/BVAMS_calibration_files/Ben_calibration_July_6_2024/bluePrimaryJuly0624_initialPositionFocus3_100.mat');
spectralAxisInterp = 380:780;
energyInterp = interp1(spectralAxis,energy,spectralAxisInterp,'spline');
figure; plot(spectralAxisInterp,energyInterp,'b-');

