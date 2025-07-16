%% GET CHROMATICITY COORDINATES OF EXP 1 CONDITIONS

% Setting up display properties
d = displayCreate('OLED-Samsung');
d = displaySet(d, 'name', 'my display');
d = displaySet(d,'ViewingDistance',1); % simulated screen distance
d = displaySet(d,'dpi',378); % simulated screen distance

if strcmp(getenv('USER'),'benjaminchin')
    calPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/BVAMS_calibration_files/Ben_calibration_July_6_2024/';
     stimPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/';
     savePath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImages/S';
end

load([calPath 'redPrimaryJuly0624_initialPositionFocus3_100.mat']);
d.spd(:,1) = energy;
load([calPath 'greenPrimaryJuly0624_initialPositionFocus3_100.mat']);
d.spd(:,2) = energy;
load([calPath 'bluePrimaryJuly0624_initialPositionFocus3_100.mat']);
d.spd(:,3) = energy;
d.gamma(:,1) = (d.gamma(:,1).^(1/2.2)).^2.5;
d.gamma(:,2) = (d.gamma(:,2).^(1/2.2)).^2.7;
d.gamma(:,3) = (d.gamma(:,3).^(1/2.2)).^2.3;

% COLOR MATCHING FUNCTIONS
S = [380 4 101]; % weird convention used by Brainard lab for defining wavelengths
load T_xyz1931; % load color matching functions
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale
wave = S(1):S(2):S(1)+S(2)*(S(3)-1); % define wavelength vector

colorsAll = [0.327 0.000 1.000; ...
             0.432 0.000 1.000; ...
             0.569 0.000 1.000; ...
             0.569 0.000 0.740; ...
             0.569 0.000 0.547; ...
             0.327 0.334 1.000; ...
             0.432 0.334 1.000; ...
             0.569 0.334 1.000; ...
             0.569 0.334 0.740; ...
             0.569 0.334 0.547; ...
             0.569 0.432 1.000; ...  
            ];

for i = 1:size(colorsAll,1)
    rgbGammaCorrected = [colorsAll(i,1)^2.5 colorsAll(i,2)^2.7 colorsAll(i,3)^2.3];
    spdTmp = d.spd(:,1).*rgbGammaCorrected(1) + d.spd(:,2).*rgbGammaCorrected(2) + d.spd(:,3).*rgbGammaCorrected(3);
    colorCndXYZ(:,i) = T_sensorXYZ*spdTmp;
    colorCndxyY(:,i) = XYZToxyY(colorCndXYZ(:,i));
    colorCndHSV(:,i) = rgb2hsv(colorsAll(i,:));
end

figure;
plotChromaticity;
hold on;
plot(colorCndxyY(1,1:5),colorCndxyY(2,1:5),'k^','LineWidth',1,'MarkerSize',10);
plot(colorCndxyY(1,1:5),colorCndxyY(2,6:10),'ks','LineWidth',1,'MarkerSize',10);
plot(colorCndxyY(1,11),colorCndxyY(2,11),'ko','LineWidth',1,'MarkerSize',10);
legend('','No green','With green','All equal lum');

%%

clear all; 
