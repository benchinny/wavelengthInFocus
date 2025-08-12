%% INITIAL PARAMETERS

clear;

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

subjNum = 20;

%% MEAN COEFFICIENTS FROM ACCOMMODATION EXP

[meanC, rgb1all, meanv00all] = ARCnlz_mainExpSortColorAbb(subjNum+10,dataPath);

indPurple = abs(rgb1all(:,1)-0.569)<0.001 & ...
            abs(rgb1all(:,2)-0.000)<0.001 & ...
            abs(rgb1all(:,3)-1.000)<0.001 & ...
            abs(meanv00all-2.5)<0.001;
meanCpurple = meanC(indPurple,:);

%% DATA FROM ACUITY EXP

wvfFiles = ARCacuAnalysisWvfSubj(subjNum, dataPath);

dataFolder = [dataPath 'data\csvFiles\SUBJ\'];

cAll = [];

for i = 1:length(wvfFiles)
    ZernikeTable = readtable([dataFolder wvfFiles{i}]);
    NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
    c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
    PARAMS = struct;
    PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
    PARAMS.PupilFitSize=mean(table2array(ZernikeTable(:,5))); 
    PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
    c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
    cAll = [cAll; c];
end

indBad = cAll(:,4)==0 | cAll(:,4)<-10;
meanCpurpleAcu = mean(cAll(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS

%%

figure;
hold on;
plot(mean(meanCpurple(:,3:end)));
plot(meanCpurpleAcu(3:end))
axis square;