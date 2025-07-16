%% LOAD MAIN EXPERIMENT FILES

subjNum = 15;

if subjNum==11
   % blockNums = 2:7;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   blockNums = 11:16;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};   
   subjName = ['S' num2str(subjNum) '-OD'];
elseif subjNum==12
   % blockNums = 2:7;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   % subjName = ['S' num2str(subjNum) '-OD'];

   blockNums = [11:15 18];
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];
elseif subjNum==13
   % blockNums = 3:8;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};
   blockNums = 12:17;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};   
   subjName = ['S' num2str(subjNum) '-OD'];   
elseif subjNum==14
   blockNums = 9:14;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   % blockNums = 3:8;
   % trialNums = {1:33 1:33 1:33 1:33 1:33 1:33};   
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==15
   blockNums = 3:8;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];   
elseif subjNum==17
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==18
   blockNums = 10:15;
   % blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==19
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];         
elseif subjNum==20
   blockNums = 3:8;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];         
elseif subjNum==21
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];            
elseif subjNum==22
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];          
elseif subjNum==26
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];      
elseif subjNum==27
   blockNums = 2:7;
   trialNums = {1:36 1:36 1:36 1:36 1:36 1:36};
   subjName = ['S' num2str(subjNum) '-OD'];         
end

meanC = [];
rgb1all = [];
meanv00all = [];
nIndBadTracker = [];
c4all = {};

for i = 1:length(blockNums)
    blockNumTmp = blockNums(i);
    trialNumsTmp = trialNums{i};
    AFCp = ARCloadFileBVAMS(subjNum,blockNumTmp);
    for j = 1:length(trialNumsTmp)
        [ZernikeTable, ~, ~, ~] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumsTmp(j),0);
        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65.
        indBadPupil = table2array(ZernikeTable(:,5))==0;
        PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(~indBadPupil,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        % THESE ARE ALL THE ZERNIKE COEFFICIENTS FOR THIS TRIAL. THERE ARE
        % 65 OF THEM. BUT ONLY THE 4TH CORRESPONDS TO DEFOCUS
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        % THE THREE LINES BELOW ARE FOR REMOVING BLINKS
        indBad = c(:,4)==0;
        nIndBadTracker(end+1) = sum(indBad);
        c(indBad,4) = mean(c(~indBad,4));
        % TAKE MEAN OF COEFFICIENTS 
        meanC(end+1,:) = mean(c,1);    
        % STORE THEM
        c4all{end+1} = c(:,4);
    end
    % RGB VALUES FOR EACH TRIAL
    rgb1all = [rgb1all; AFCp.rgb100];
    % STIMULUS OPTICAL DISTANCE FOR EACH TRIAL
    meanv00all = [meanv00all; AFCp.meanv00./1.2255];
end

%% PLOTTING ALL TRIALS FOCUS FOR EACH OF 3 OPTICAL DISTANCES

% MAKE SURE YOU DIVIDE BY THIS VALUE TO CONVERT THE DEFOCUS TERM FROM THE
% WAVEFRONT SENSOR TO THE ACTUAL DEFOCUS VALUE
defocusScaleFactor = 1./((1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2));
defocusAddFactor = 0.9274;

% DEFOCUS AT 550NM. REMEMBER THE WAVEFRONT SENSOR MEASURES THE WAVEFRONT AT
% 875NM, SO TO PREDICT THE DEFOCUS AT A REASONABLE WAVELENGTH SUCH AS
% 550NM, WE NEED TO CORRECT BY THE EXPECTED AMOUNT OF LONGITUDINAL
% CHROMATIC ABERRATION
defocusAt550 = defocusAddFactor+meanC(:,4).*defocusScaleFactor;

figure; 
plot(meanv00all,defocusAt550,'ko');
set(gca,'FontSize',15);
xlim([1 4]);
axis square;
xlabel('Stimulus optical distance');
ylabel('Raw refractive power (D)');

%% PLOTTING ALL TRIAL MEANS PER CONDITION AND DISTANCE

% THESE TWO PARAMETERS ARE REQUIRED FOR CONVERTING RGB VALUES TO LUMINANCE
% VALUES
lumScaleRGB = [4.0888 9.6669 1];
gammaRGB = [2.5 2.7 2.3];

% LUMINANCE OF THE RED, GREEN, AND BLUE PRIMARIES OF THE DISPLAY
rgbLumNorm = [lumScaleRGB(1).*rgb1all(:,1).^gammaRGB(1) lumScaleRGB(2).*rgb1all(:,2).^gammaRGB(2) lumScaleRGB(3).*rgb1all(:,3).^gammaRGB(3)];

% PRE-DEFINING COLOR CONDITIONS TO LOOK FOR
conditionsOrderedNormRGB = [0.25 0.00 1.00; ...
                            0.50 0.00 1.00; ...
                            1.00 0.00 1.00; ...
                            1.00 0.00 0.50; ...
                            1.00 0.00 0.25; ...
                            0.25 0.50 1.00; ...
                            0.50 0.50 1.00; ...
                            1.00 0.50 1.00; ...
                            1.00 0.50 0.50; ...
                            1.00 0.50 0.25];

% OPTICAL DISTANCES TO ANALYZE
optDistToCheckAll = [1.5 2.5 3.5];

%%

% UNIQUE RED GREEN BLUE VALUES (IN TERMS OF NORMALIZED LUMINANCE)
rgbLumNormUnq = unique(rgbLumNorm,'rows');
% INDEX TO GET ALL TRIALS NUMBERS WITH SAME COLOR AND OPTICAL DISTANCE
ind = find(rgbLumNorm(:,1)==rgbLumNormUnq(1,1) & rgbLumNorm(:,2)==rgbLumNormUnq(1,2) & rgbLumNorm(:,3)==rgbLumNormUnq(1,3) & meanv00all==1.5);

cTmp = []; % INITIALIZING VECTOR CONTAINING ALL MEASUREMENTS FOR A GIVEN CONDITION
for j = 1:length(ind) % LOOPING OVER TRIAL INDICES
    cTmp = [cTmp; c4all{ind(j)}]; % COMBINING ALL TRIALS IN THAT CONDITION INTO ONE VECTOR
end

% std(cTmp.*defocusScaleFactor+defocusAddFactor)
% figure; hist(cTmp.*defocusScaleFactor+defocusAddFactor,21)
