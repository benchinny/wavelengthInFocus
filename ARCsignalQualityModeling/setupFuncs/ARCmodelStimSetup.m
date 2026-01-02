function [s1, s2] = ARCmodelStimSetup(dataPath,subjNum,accommodationORacuity,d,rgb,bPLOT)

% makes stimulus struct for retinal image modeling (for both accommodation
% and acuity experiments)
%
% dataPath: path to all data
% subjNum: subject number
% accommodationORacuity: which experiment we are modeling
%                        'accommodation' -> experiment 1
%                        'acuity'-> experiment 2
% d: display struct (must be made in ISETBIO using ARCmodelDispSetup)
% rgb: stimulus color
% bPLOT: plot or not
%
% s1: output stimulus struct. If modeling accommodation experiment, this is
%     a three-letter word. If modeling acuity experiment, this is a Gabor
%     oriented 15 degrees. 
% s2: output stimulus struct. If modeling accommodation experiment, this is
%     empty. If modeling acuity experiment, this is a Gabor oriented -15 degrees. 

% GAMMA (BASED ON BVAMS CALIBRATION)
gammaR = 2.5;
gammaG = 2.7;
gammaB = 2.3;

if strcmp(accommodationORacuity,'accommodation')
    % PATH TO STIMULUS SPATIAL PATTERN
    stimPath = fullfile(dataPath,'data','helperFiles');
    % READ IN STIMULUS SPATIAL PATTERN
    im = imread(fullfile(stimPath,'word_image_01.png'));
    im = double(im); % CONVERT TO DOUBLE FOR INCREASED PRECISION)
    % IN THE IMAGE THE BLUE CHANNEL IS AT MAXIMUM VALUE, SO WE USE THIS
    % PATTERN AS THE 'TEMPLATE'
    imPattern = squeeze(im(:,:,3));
    % ADD ZERO PADDING (TO AVOID ARTIFACTS IN LATER ANALYSES)
    imPattern = [zeros([100 size(imPattern,2)]); imPattern; zeros([100 size(imPattern,2)])];
    imPattern = [zeros([size(imPattern,1) 30]) imPattern zeros([size(imPattern,1) 30])];
    % NOW FILL IN ALL CHANNELS WITH SPATIAL PATTERN SCALED BY RGB VALUES
    I(:,:,3) = rgb(3).*imPattern;
    I(:,:,2) = rgb(2).*imPattern;
    I(:,:,1) = rgb(1).*imPattern;
    % NORMALIZE TO 1 (WHAT ISETBIO EXPECTS)
    I = I./255;
    
    % Turn image into 'scene'
    s1 = sceneFromFile(I, 'rgb', [], d);  % The display is included here
    % I think this copies the struct into an object
    vcAddObject(s1); 
    % FOR ACCOMMODATION EXPERIMENT, THERE IS NO SECOND STIMULUS
    s2 = [];
    
    % SOMETIMES MIGHT WANT TO PLOT THE STIMULUS TO MAKE SURE NOTHING GOT
    % MESSED UP DURING THE STIMULUS CREATION PROCESS
    if bPLOT 
        figure; 
        set(gcf,'Position',[289 428 1056 420]);
        subplot(1,2,2);
        imagesc(I);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        axis square;
        set(gca,'FontSize',15);
        title('Original');
        subplot(1,2,2);
        plot(s1.spectrum.wave,squeeze(s1.data.photons(160,160,:)),'-k','LineWidth',1);
        formatFigure('Wavelength (\lambda)','Photons');
        axis square;
    end
elseif strcmp(accommodationORacuity,'acuity')
    % GET THRESHOLDS FROM CONTRAST CALIBRATION EXPERIMENT
    thresholds = ARCnlz_contrastThresholds(subjNum,0,dataPath);    
    % GABOR PARAMETERS
    frqCpd = 15;
    contrast = thresholds(4);
    stimPositionsX = smpPos(260,390); % SAMPLING FOR MAKING GABOR
    x0 = 0; % CENTER OF GABOR
    y0 = 0; % CENTER OF GABOR
    frqCpdAll = [frqCpd 3*frqCpd 5*frqCpd 7*frqCpd]; % ALL FREQUENCIES FOR SQUARE WAVE
    contrastAll = [contrast contrast/3 contrast/5 contrast/7]; % ALL CONTRASTS FOR SQUARE WAVE
    orientation1 = 15; % STIMULUS 1 ORIENTATION
    orientation2 = -15; % STIMULUS 2 ORIENTATION
    phs = 90; % PHASE
    sigmaX = 0.2; % STANDARD DEVIATION OF GAUSSIAN ENVELOPE IN X
    sigmaY = 0.2; % STANDARD DEVIATION OF GAUSSIAN ENVELOPE IN Y
    
    % MODEL ACUITY STIMULUS FOR +15 DEG ORIENTATION
    acuStimOrig1 = ARC2Dgabor(stimPositionsX,[],x0,y0,frqCpdAll, ...
                   contrastAll,orientation1,phs,sigmaX,sigmaY, ...
                   [rgb(1,1)^gammaR rgb(1,2)^gammaG rgb(1,3)^gammaB],1,1,0,0);
    
    % MAKE SURE TO GAMMA CORRECT STIMULUS FOR MODELING
    acuStimOrig1(:,:,1) = acuStimOrig1(:,:,1).^(1/gammaR);
    acuStimOrig1(:,:,2) = acuStimOrig1(:,:,2).^(1/gammaG);
    acuStimOrig1(:,:,3) = acuStimOrig1(:,:,3).^(1/gammaB);
    I1 = acuStimOrig1;
    
    % MODEL ACUITY STIMULUS FOR -15 DEG ORIENTATION
    acuStimOrig2 = ARC2Dgabor(stimPositionsX,[],x0,y0,frqCpdAll, ...
                   contrastAll,orientation2,phs,sigmaX,sigmaY, ...
                   [rgb(1,1)^gammaR rgb(1,2)^gammaG rgb(1,3)^gammaB],1,1,0,0);
    
    % MAKE SURE TO GAMMA CORRECT STIMULUS FOR MODELING
    acuStimOrig2(:,:,1) = acuStimOrig2(:,:,1).^(1/gammaR);
    acuStimOrig2(:,:,2) = acuStimOrig2(:,:,2).^(1/gammaG);
    acuStimOrig2(:,:,3) = acuStimOrig2(:,:,3).^(1/gammaB);
    I2 = acuStimOrig2;
    
    % Turn image into 'scene'
    s1 = sceneFromFile(I1, 'rgb', [], d);  % The display is included here
    s2 = sceneFromFile(I2, 'rgb', [], d);
    % I think this copies the struct into an object
    vcAddObject(s1); 
    vcAddObject(s2); 
else
    error('ARCmodelStimSetup: invalid value of accommodationORacuity input parameter. This function only handles the accommodation or acuity experiments!');
end

end