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

bUseBVAMScal = 1; % if using BVAMS calibration data

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

if bUseBVAMScal % LOAD CALIBRATION DATA
    load([calPath 'redPrimaryJuly0624_initialPositionFocus3_100.mat']);
    d.spd(:,1) = energy;
    load([calPath 'greenPrimaryJuly0624_initialPositionFocus3_100.mat']);
    d.spd(:,2) = energy;
    load([calPath 'bluePrimaryJuly0624_initialPositionFocus3_100.mat']);
    d.spd(:,3) = energy;
end
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

%%

if subjNum==3
   blockNums = 12:17;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+10) '-OD'];
elseif subjNum==10
   blockNums = 3:8;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
   subjName = ['S' num2str(subjNum+10) '-OD'];
elseif subjNum==1
   blockNums = 11:16;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
   subjName = ['S' num2str(subjNum+10) '-OD'];
elseif subjNum==5
   blockNums = 3:8;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
   subjName = ['S' num2str(subjNum+10) '-OD'];
elseif subjNum==9
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
   subjName = ['S' num2str(subjNum+10) '-OD'];
elseif subjNum==16
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+10) '-OD'];
elseif subjNum==17
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+10) '-OD'];
elseif subjNum==18
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+10) '-OD'];
elseif subjNum==20
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+10) '-OD'];
end

%%

cAll = [];
optDistAll = [];
rgbAll = [];

for l = 1:length(blockNums) % LOOP OVER BLOCK
    blockNumInd = l;
    blockNumTmp = blockNums(blockNumInd);
    AFCp = ARCloadFileBVAMS(subjNum+10,blockNumTmp,dataPath); % LOAD BVAMS DATA
    rgbAll = [rgbAll; AFCp.rgb100];
    for k = 1:36 % LOOP OVER TRIAL
        trialNumTmp = trialNums(k,blockNumInd);
        
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumTmp,0,dataPath);

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS = struct;
        PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilSize = 4.*ones(size(PARAMS.PupilSize)); % hard code pupil size (FOR BVAMS ONLY)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(:,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        optDistTmp = (AFCp.meanv00(k)./1.2255).*ones([size(c,1) 1]);
        optDistAll = [optDistAll; optDistTmp];
        cAll = [cAll; c];
    end
end

indBad = cAll(:,4)==0;
meanC = mean(cAll(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS
rgb00 = unique(rgbAll,'rows');

for k = 1 % LOOP OVER TRIAL
    % recreate stimulus
    rVal = rgb00(k,1);
    gVal = rgb00(k,2);
    bVal = rgb00(k,3);
    im = imread([stimPath '/word_image_01.png']);
    im = double(im);
    imPattern = squeeze(im(:,:,3));
    imPattern = [zeros([100 size(imPattern,2)]); imPattern; zeros([100 size(imPattern,2)])];
    imPattern = [zeros([size(imPattern,1) 30]) imPattern zeros([size(imPattern,1) 30])];
    I(:,:,3) = bVal.*imPattern;
    I(:,:,2) = gVal.*imPattern;
    I(:,:,1) = rVal.*imPattern;
    I = I./255;
    
    % Turn image into 'scene'
    s = sceneFromFile(I, 'rgb', [], d);  % The display is included here
    % I think this copies the struct into an object
    vcAddObject(s); 
    
    % figure; 
    % set(gcf,'Position',[289 428 1056 420]);
    % subplot(1,3,1);
    % plot(d.wave,d.spd(:,1),'r','LineWidth',1.5); hold on;
    % plot(d.wave,d.spd(:,2),'g','LineWidth',1.5);
    % plot(d.wave,d.spd(:,3),'b','LineWidth',1.5);
    % axis square;
    % formatFigure('Wavelength (\lambda)','Radiance');
    % subplot(1,3,2);
    % imagesc(I);
    % set(gca,'XTick',[]);
    % set(gca,'YTick',[]);
    % axis square;
    % set(gca,'FontSize',15);
    % title('Original');
    % subplot(1,3,3);
    % plot(s.spectrum.wave,squeeze(s.data.photons(160,160,:)),'-k','LineWidth',1);
    % formatFigure('Wavelength (\lambda)','Photons');
    % axis square;
    
    % Computing peak correlation for different wavelengths in focus
    
    wave2 = 380:4:780;

    for i = 41
        % zCoeffs = [0 zeros(size(meanC(1:end-1)))];
        zCoeffs = [0 meanC(1:end-1)];
        wvfP = wvfCreate('calc wavelengths', wave, ...
            'measured wavelength', wave2(i), ...
            'zcoeffs', zCoeffs, 'measured pupil', PARAMS.PupilSize, ...
            'name', sprintf('human-%d', PARAMS.PupilSize),'spatial samples',size(I,2));
        wvfP.calcpupilMM = PARAMS.PupilSize;
        if subjNum==1
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
        if defocusFromLCA<1
            wvfP.refSizeOfFieldMM = 12;
        else
            wvfP.refSizeOfFieldMM = 6;
        end
        wvfP = wvfSet(wvfP, 'zcoeff', 0, 'defocus');
        
        % Convert to siData format as well as wavefront object
        [siPSFData, wvfP] = wvf2SiPsfARC(wvfP,'showBar',false,'nPSFSamples',size(I,2),'umPerSample',1.1512); 
        oi = wvf2oi(wvfP); % CONVERT TO OPTICS OBJECT
        % PADDING PSF NECESSARY TO ENSURE SAME SIZE AS STIMULUS
        paddingXCpsf = round((size(siPSFData.psf,2)-size(s.data.photons,2))/2);
        paddingYRpsf = round((size(siPSFData.psf,1)-size(s.data.photons,1))/2);
        indNotPadded = {(paddingYRpsf+1):(size(siPSFData.psf,1)-paddingYRpsf) ...
                        (paddingXCpsf+1):(size(siPSFData.psf,2)-paddingXCpsf)};
        % I HAD TO WRITE NEW CODE TO SET PSF BECAUSE I COULDN'T FIGURE OUT
        % HOW TO DO WHAT I WANTED WITHIN ISETBIO FRAMEWORK
        oi.optics.OTF = [];
        for j = 1:size(siPSFData.psf,3)
            oi.optics.OTF.OTF(:,:,j) = fft2(fftshift(squeeze(siPSFData.psf(indNotPadded{1},indNotPadded{2},j))));
        end
        % oi = oiCreateARC('human',wave,Dall(i)); % create optics
        oi = oiCompute(oi, s); % compute optical image of stimulus
    
        % Create the coneMosaic object
        cMosaic = coneMosaic;

        % Set size to show relevant portion of scene
        cMosaic.setSizeToFOV(1 * sceneGet(s, 'fov'));

        % key line for computing absorptions
        absorptions = cMosaic.computeSingleFrame(oi, 'fullLMS', true);            
        
        display(['Peak correlation loop ' num2str(i) ' stimulus ' num2str(k)]);
        absorptions = single(absorptions);
        S = struct;
        S.absorptions = absorptions;
        fnameCone = ['subj' num2str(subjNum) 'stimulus' num2str(k) 'focusInd' num2str(i)];
        % save([savePath num2str(subjNum) '/' fnameCone '.mat'],"-fromstruct",S);
    end
end

end
% save('/Users/benjaminchin/Documents/ARchromaScraps/ARCmodelOutput2_5.mat','meanv00all','wvInFocus1all','rgb1all','defocusBasic');
