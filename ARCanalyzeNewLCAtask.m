%%

% subj = 'BenChin-OD';
subj = 'S12-OD';
if strcmp(subj,'BenChin-OD')
    vs = 12;
end
if strcmp(subj,'S12-OD')
    vs = 17;
end

[ZernikeTable, ~, ~, ~, trialTagIndex] = ARCloadFileFIATallInstances(subj,vs,0);
NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65.
indBadPupil = table2array(ZernikeTable(:,5))==0;
PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
PARAMS.PupilFitSize=mean(table2array(ZernikeTable(~indBadPupil,5))); 
PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));

indBad = c(:,4)==0;
cGood = c(~indBad,4);
defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
defocusAt550 = humanWaveDefocus(875)-humanWaveDefocus(550)+cGood./defocusCorrectionFactor;
trialTagIndexGood = trialTagIndex(~indBad);

% LOAD LCA FILES

if strcmp(subj,'BenChin-OD')
    load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/ARC/LCAfileS1.mat');
    clrIndAll = clrIndAll(1:30)';
    powerDispRall = powerDispRall(1:30);
end
if strcmp(subj,'S12-OD')
    load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/ARC/LCAadjustmentS12.mat');
end

for i = 1:length(unique(trialTagIndexGood))
    trialInd = trialTagIndexGood==i;
    defocusAt550mean(i) = mean(defocusAt550(trialInd));
end

mean(0.816*(14.3-powerDispRall(clrIndAll'==1 & powerDispRall'<15)')-defocusAt550mean(clrIndAll'==1 & powerDispRall'<15)')
mean(0.816*(14.3-powerDispRall(clrIndAll'==2 & powerDispRall'<15)')-defocusAt550mean(clrIndAll'==2 & powerDispRall'<15)')
mean(0.816*(14.3-powerDispRall(clrIndAll'==3 & powerDispRall'<15)')-defocusAt550mean(clrIndAll'==3 & powerDispRall'<15)')
