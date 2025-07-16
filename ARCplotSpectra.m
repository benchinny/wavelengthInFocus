%%

figure;
hold on;
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/BVAMS_calibration_files/Ben_calibration_July_6_2024/redPrimaryJuly0624_initialPositionFocus3_100.mat');
plot(spectralAxis,energy./max(energy),'-r','LineWidth',1);
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/BVAMS_calibration_files/Ben_calibration_July_6_2024/bluePrimaryJuly0624_initialPositionFocus3_100.mat');
plot(spectralAxis,energy./max(energy),'-b','LineWidth',1);
load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/BVAMS_calibration_files/Ben_calibration_July_6_2024/greenPrimaryJuly0624_initialPositionFocus3_100.mat');
plot(spectralAxis,energy./max(energy),'-g','LineWidth',1);
axis square;
set(gca,'FontSize',15);
xlabel('Wavelength (nm)');
ylabel('Normalized Power');

%%

clear all;