function oi = ARCimgQualityAnalysisConesMeanZ(subjNum,dataPath)

% subjNum values for participants who passed screening: 1, 3, 5, 10, 16,
% 17, 18, 20

%% Initialize and clear
ieInit;

%% Set up display struct and build Ben's stimulus

% Setting up display properties
d = displayCreate('OLED-Samsung');
d = displaySet(d, 'name', 'my display');
% SIMULATED SCREEN DISTANCE--IT'S ARBITRARY SINCE THE BVAMS HAS
% NON-STANDARD OPTICAL PROPERTIES
d = displaySet(d,'ViewingDistance',1);
% 378 DOTS PER INCH YIELDS APPROXIMATELY 260 PIXELS PER VISUAL DEGREE,
% WHICH IS WHAT THE BVAMS HAS
d = displaySet(d,'dpi',378); 
% GET RID OF ALL UNNECESSARY FIELDS
d.dixel = [];
d.mainimage = [];

bPlotStim = false;

if ispc
    slash = '\';
else
    slash = '/';
end

% PATH TO CALIBRATION DATA
calPath = [dataPath 'BVAMS_calibration_files' slash 'Ben_calibration_July_6_2024' slash];
% PATH TO STIMULUS SPATIAL PATTERN
stimPath = [dataPath 'stimuli' slash];
% PATH TO SAVE
savePath = [dataPath 'data' slash 'coneImages' slash 'S'];

% LOAD BVAMS CALIBRATION DATA
load([calPath 'redPrimaryJuly0624_initialPositionFocus3_100.mat']);
d.spd(:,1) = energy;
load([calPath 'greenPrimaryJuly0624_initialPositionFocus3_100.mat']);
d.spd(:,2) = energy;
load([calPath 'bluePrimaryJuly0624_initialPositionFocus3_100.mat']);
d.spd(:,3) = energy;

% ISETBIO DISPLAY STRUCT HAS DEFAULT GAMMA OF 2.2, SO NEED TO UNDO IT, THEN
% APPLY OUR EMPIRICALLY DERIVED GAMMA FROM CALIBRATION MEASUREMENTS
d.gamma(:,1) = (d.gamma(:,1).^(1/2.2)).^2.5;
d.gamma(:,2) = (d.gamma(:,2).^(1/2.2)).^2.7;
d.gamma(:,3) = (d.gamma(:,3).^(1/2.2)).^2.3;

% COLOR MATCHING FUNCTIONS
S = [380 4 101]; % weird convention used by Brainard lab for defining wavelengths
load T_xyz1931; % load color matching functions
% DEFINE WAVELENGTH VECTOR: 380 TO 780 WITH 4NM INCREMENTS
wave = S(1):S(2):S(1)+S(2)*(S(3)-1);

% FILE NUMBERS FOR DIFFERENT SUBJECTS (CORRESPONDING TO EACH BLOCK)
if subjNum==3
   blockNums = 12:17;
elseif subjNum==10
   blockNums = 3:8;
elseif subjNum==1
   blockNums = 11:16;
elseif subjNum==5
   blockNums = 3:8;
elseif subjNum==9
   blockNums = 2:7;
elseif subjNum==16
   blockNums = 2:7;
elseif subjNum==17
   blockNums = 2:7;
elseif subjNum==18
   blockNums = 2:7;
elseif subjNum==20
   blockNums = 2:7;
end

% ALL SUBJECTS HAD THE SAME NUMBER OF TRIALS
trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
% FOR LOADING DATA LATER
subjName = ['S' num2str(subjNum+10) '-OD'];

%%

cAll = []; % INITIALIZE MATRIX FOR STORING COEFFICIENTS
optDistAll = []; % FOR CONCATENATING OPTICAL DISTANCES PER TRIAL
rgbAll = []; % FOR CONCATENATING RGB VALUES PER TRIAL

% GET THE AVERAGE HIGHER-ORDER ABERRATIONS FOR THIS SUBJECT ACROSS ALL
% TRIALS
for l = 1:length(blockNums) % LOOP OVER BLOCK
    blockNumInd = l; % HELPS WITH READABILITY
    blockNumTmp = blockNums(blockNumInd); % GRAB BLOCK
    AFCp = ARCloadFileBVAMS(subjNum+10,blockNumTmp,dataPath); % LOAD BVAMS DATA
    rgbAll = [rgbAll; AFCp.rgb100]; % STACK COLOR CONDITIONS
    for k = 1:36 % LOOP OVER TRIAL
        trialNumTmp = trialNums(k,blockNumInd); % GRAB TRIAL NUMBERS
        
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumTmp,0,dataPath);

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS = struct;
        indBadPupil = table2array(ZernikeTable(:,5))==0; % IDENTIFY BLINKS
        % STORE PUPIL SIZE (SHOULD BE 4MM FOR EVERYTHING)
        PARAMS.PupilSize=mean(table2array(ZernikeTable(~indBadPupil,5))); 
        % GRAB ALL ZERNIKE COEFFICIENTS
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        % STACK OPTICAL DISTANCES. DIVIDE BY 1.2255 TO ACCOUNT FOR BVAMS
        % CALIBRATION RELATING OPTOTUNE DEFOCUS CHANGE TO DEFOCUS AT EYE
        optDistTmp = (AFCp.meanv00(k)./1.2255).*ones([size(c,1) 1]);
        optDistAll = [optDistAll; optDistTmp];
        cAll = [cAll; c]; % STORE ALL COEFFICIENTS
    end
end

indBad = cAll(:,4)==0; % LOOK FOR BLINKS IN DEFOCUS ABERRATION VECTOR
meanC = mean(cAll(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS
rgb00 = unique(rgbAll,'rows'); % GET UNIQUE COLOR CONDITIONS

for k = 1:size(rgb00,1) % LOOP OVER COLOR CONDITIONS
    % GET RGB VALUES FOR A PARTICULAR CONDITION
    rVal = rgb00(k,1);
    gVal = rgb00(k,2);
    bVal = rgb00(k,3);
    % READ IN STIMULUS SPATIAL PATTERN
    im = imread([stimPath '/word_image_01.png']);
    im = double(im); % CONVERT TO DOUBLE FOR INCREASED PRECISION)
    % IN THE IMAGE THE BLUE CHANNEL IS AT MAXIMUM VALUE, SO WE USE THIS
    % PATTERN AS THE 'TEMPLATE'
    imPattern = squeeze(im(:,:,3));
    % ADD ZERO PADDING (TO AVOID ARTIFACTS IN LATER ANALYSES)
    imPattern = [zeros([100 size(imPattern,2)]); imPattern; zeros([100 size(imPattern,2)])];
    imPattern = [zeros([size(imPattern,1) 30]) imPattern zeros([size(imPattern,1) 30])];
    % NOW FILL IN ALL CHANNELS WITH SPATIAL PATTERN SCALED BY RGB VALUES
    I(:,:,3) = bVal.*imPattern;
    I(:,:,2) = gVal.*imPattern;
    I(:,:,1) = rVal.*imPattern;
    % NORMALIZE TO 1 (WHAT ISETBIO EXPECTS)
    I = I./255;
    
    % Turn image into 'scene'
    s = sceneFromFile(I, 'rgb', [], d);  % The display is included here
    % I think this copies the struct into an object
    vcAddObject(s); 
    
    % SOMETIMES MIGHT WANT TO PLOT THE STIMULUS TO MAKE SURE NOTHING GOT
    % MESSED UP DURING THE STIMULUS CREATION PROCESS
    if bPlotStim 
        figure; 
        set(gcf,'Position',[289 428 1056 420]);
        subplot(1,3,1);
        plot(d.wave,d.spd(:,1),'r','LineWidth',1.5); hold on;
        plot(d.wave,d.spd(:,2),'g','LineWidth',1.5);
        plot(d.wave,d.spd(:,3),'b','LineWidth',1.5);
        axis square;
        formatFigure('Wavelength (\lambda)','Radiance');
        subplot(1,3,2);
        imagesc(I);
        set(gca,'XTick',[]);
        set(gca,'YTick',[]);
        axis square;
        set(gca,'FontSize',15);
        title('Original');
        subplot(1,3,3);
        plot(s.spectrum.wave,squeeze(s.data.photons(160,160,:)),'-k','LineWidth',1);
        formatFigure('Wavelength (\lambda)','Photons');
        axis square;
    end
    
    % PEAK CORRELATION WILL BE COMPUTED FOR A RANGE OF 'WAVELENGTHS IN
    % FOCUS'. THIS IS A DISTINCT VARIABLE FROM 'wave', WHICH IS THE
    % WAVELENGTHS OVER WHICH THE STIMULUS IS DEFINED. I HAVE MADE THEM THE
    % SAME, BUT THEY DON'T NECESSARILY HAVE TO BE (E.G. IF YOU JUST WANT TO
    % LOOK WHAT HAPPENS WHEN A SPECIFIC WAVELENGTH IS IN FOCUS)
    wave2 = 380:4:780;

    for i = 1:length(wave2) % LOOP OVER WAVELENGTHS IN FOCUS
        % % IF YOU WANT TO COMPUTE FOR DIFFRACTION-LIMITED SCENARIO 
        % zCoeffs = [0 zeros(size(meanC(1:end-1)))];

        % REFORMAT COEFFICIENTS (WAVEFRONT SENSOR LEAVES OUT THE 'PISTON'
        % ZERNIKE TERM, WHICH IS ALWAYS 0)
        zCoeffs = [0 meanC(1:end-1)];
        % CREATE ISETBIO WAVEFRONT OBJECT
        wvfP = wvfCreate('calc wavelengths', wave, ...
            'measured wavelength', wave2(i), ...
            'zcoeffs', zCoeffs, 'measured pupil', PARAMS.PupilSize, ...
            'name', sprintf('human-%d', PARAMS.PupilSize),'spatial samples',size(I,2));
        % MAKE SURE THE calcpupilMM FIELD MATCHES THE ACTUAL PUPIL SIZE IN
        % THE EXPERIMENT
        wvfP.calcpupilMM = PARAMS.PupilSize;
        % SET CUSTOM LCA FUNCTION PER SUBJECT--ISETBIO WANTS IT TO BE SET
        % IN A PARTICULAR FORMAT
        if subjNum==1
            % CALCULATE MINIMUM AND MAXIMUM DEFOCUS FROM LCA FOR A
            % PARTICULAR SUBJECT FOR THE RANGE OF WAVELENGTHS TO ANALYZE.
            % THIS WILL ENSURE THAT THE MESH OVER THE PUPIL FUNCTION SPANS
            % A SUFFICIENT RANGE (refSizeOfFieldMM). IF DEFOCUS IS LARGE
            % ENOUGH, THE RANGE NEEDS TO BE REDUCED.
            defocusFromLCA = max(abs([humanWaveDefocusS1(wave2(i),min(wave)) ...
                                      humanWaveDefocusS1(wave2(i),max(wave))]));
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS1);
        elseif subjNum==3
            defocusFromLCA = max(abs([humanWaveDefocusS3(wave2(i),min(wave)) ...
                                      humanWaveDefocusS3(wave2(i),max(wave))]));  
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS3);
        elseif subjNum==5
            defocusFromLCA = max(abs([humanWaveDefocusS5(wave2(i),min(wave)) ...
                                      humanWaveDefocusS5(wave2(i),max(wave))]));  
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS5); 
        elseif subjNum==9
            defocusFromLCA = max(abs([humanWaveDefocusS9(wave2(i),min(wave)) ...
                                      humanWaveDefocusS9(wave2(i),max(wave))]));  
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS9); 
        elseif subjNum==10
            defocusFromLCA = max(abs([humanWaveDefocusS10(wave2(i),min(wave)) ...
                                      humanWaveDefocusS10(wave2(i),max(wave))]));  
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS10); 
        elseif subjNum==16
            defocusFromLCA = max(abs([humanWaveDefocusS16(wave2(i),min(wave)) ...
                                      humanWaveDefocusS16(wave2(i),max(wave))]));  
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS16); 
        elseif subjNum==17
            defocusFromLCA = max(abs([humanWaveDefocusS17(wave2(i),min(wave)) ...
                                      humanWaveDefocusS17(wave2(i),max(wave))]));  
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS17); 
        elseif subjNum==18
            defocusFromLCA = max(abs([humanWaveDefocusS18(wave2(i),min(wave)) ...
                                      humanWaveDefocusS18(wave2(i),max(wave))]));  
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS18); 
        elseif subjNum==20
            defocusFromLCA = max(abs([humanWaveDefocusS20(wave2(i),min(wave)) ...
                                      humanWaveDefocusS20(wave2(i),max(wave))]));  
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS20); 
        else
            error('Subject number with no LCA function?');
        end

        % IF DEFOCUS IS LARGE ENOUGH, THE AREA OF THE PUPIL THE PSF IS
        % CALCULATED FROM NEEDS TO BE REDUCED, OR YOU WILL END UP WITH A
        % DEGENERATE PSF
        if defocusFromLCA<1
            wvfP.refSizeOfFieldMM = 12;
        else
            wvfP.refSizeOfFieldMM = 6;
        end
        % SET THE COEFFICIENT ON DEFOCUS TO 0 SINCE WE ARE LOOPING OVER
        % WAVELENGTHS IN FOCUS EXPLICITLY, NOT IMPLICITLY THROUGH THE
        % DEFOCUS TERM
        wvfP = wvfSet(wvfP, 'zcoeff', 0, 'defocus');
        
        % MAKE POINT-SPREAD FUNCTION (siPSFData) AND WAVEFRONT STRUCT
        [siPSFData, wvfP] = wvf2SiPsfARC(wvfP,'showBar',false,'nPSFSamples',size(I,2),'umPerSample',1.1512); 
        oi = wvf2oi(wvfP); % CONVERT WAVEFRONT STRUCT TO OPTICS OBJECT
        % NEED TO REMOVE PADDED ZEROS FROM PSF TO MAKE SAME SIZE AS
        % STIMULUS IMAGE. THE LINES BELOW IDENTIFY THE 'GOOD INDICES', I.E.
        % THE INDICES THAT AREN'T THE PADDED ZEROS TO BE REMOVED
        paddingXCpsf = round((size(siPSFData.psf,2)-size(s.data.photons,2))/2);
        paddingYRpsf = round((size(siPSFData.psf,1)-size(s.data.photons,1))/2);
        indNotPadded = {(paddingYRpsf+1):(size(siPSFData.psf,1)-paddingYRpsf) ...
                        (paddingXCpsf+1):(size(siPSFData.psf,2)-paddingXCpsf)};
        % REPLACE OTF FIELD WITH A NEW OTF CALCULATED FROM THE PSFS WE
        % JUST CALCULATED. COULDN'T FIGURE OUT HOW TO MAKE ISETBIO DO
        % THIS AUTOMATICALLY, SO I DID IT MANUALLY (AND CHECKED THE
        % OUTPUT!)
        oi.optics.OTF = []; % INITIALIZE ARRAY FOR STORING OTFS
        for j = 1:size(siPSFData.psf,3) % LOOP OVER WAVELENGTHS
            oi.optics.OTF.OTF(:,:,j) = fft2(fftshift(squeeze(siPSFData.psf(indNotPadded{1},indNotPadded{2},j))));
        end

        % MAIN STEP: COMPUTE OPTICAL IMAGE OF STIMULUS
        oi = oiCompute(oi, s); 
    
        % Create the coneMosaic object
        cMosaic = coneMosaic;

        % Set size to show relevant portion of scene
        cMosaic.setSizeToFOV(1 * sceneGet(s, 'fov'));

        % key line for computing absorptions
        absorptions = cMosaic.computeSingleFrame(oi, 'fullLMS', true);            
        
        display(['Peak correlation loop ' num2str(i) ' stimulus ' num2str(k)]);
        % CONVERT TO SINGLE TO SAVE SPACE
        absorptions = single(absorptions);
        % PLACE ABSORPTIONS IN STRUCT
        S = struct;
        S.absorptions = absorptions;
        % NAME FOR SAVING
        fnameCone = ['subj' num2str(subjNum) 'stimulus' num2str(k) 'focusInd' num2str(i)];
        % UNCOMMENT LINE BELOW TO SAVE NEW CONE IMAGES
        % save([savePath num2str(subjNum) '/' fnameCone '.mat'],"-fromstruct",S);
    end
end

end
