%%

clear; 

%%

S = [380 4 101]; % weird convention used by Brainard lab for defining wavelengths
wavelengthAxis = SToWls(S);
load('T_cones_ss2.mat');
coneFundamentals = SplineCmf(S_cones_ss2, T_cones_ss2, WlsToS(wavelengthAxis));
load T_xyz1931; % load color matching functions
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale

%%

load('/Users/benjaminchin/Documents/ARchromaScraps/testConeFundamentals.mat');
coneFundamentals(isnan(coneFundamentals))=0;

coef = [coneFundamentals(:,1:2)]\T_sensorXYZ(2,:)';

% coef = lsqlin(coneFundamentals(1:2,:)',T_sensorXYZ(2,:)',[],[],[],[],[0; 0],[1000; 1000]);

figure;
hold on;
plot(wavelengthAxis,coef(1).*coneFundamentals(:,1)+coef(2).*coneFundamentals(:,2));
plot(wavelengthAxis,T_sensorXYZ(2,:));

