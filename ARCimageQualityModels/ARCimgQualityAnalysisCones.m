%% Initialize and clear
ieInit;

%% Set up display struct and build Ben's stimulus

subjNum = 10;
subjNumEncode = subjNum+10;

% Setting up display properties
d = displayCreate('OLED-Samsung');
d = displaySet(d, 'name', 'my display');
d = displaySet(d,'ViewingDistance',1); % simulated screen distance
d = displaySet(d,'dpi',378); % simulated screen distance

bUseBVAMScal = 1; % if using BVAMS calibration data

if strcmp(getenv('USER'),'benjaminchin')
    calPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/BVAMS_calibration_files/display calibration on August3/';
     stimPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/';
     savePath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImages/S';
end

if strcmp(getenv('USER'),'benchin')
    calPath = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/BVAMS_calibration_files/display calibration on August3/';
    stimPath = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/stimuli/';
    savePath = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImages/S';
end

if bUseBVAMScal
    load([calPath 'Right_disp_Red.mat']);
    d.spd(:,1) = CurrentSpectrum.Spectral.emission_data;
    load([calPath 'Right_disp_Green.mat']);
    d.spd(:,2) = CurrentSpectrum.Spectral.emission_data;
    load([calPath 'Right_disp_Blue.mat']);
    d.spd(:,3) = CurrentSpectrum.Spectral.emission_data;
end
d.gamma(:,1) = (d.gamma(:,1).^(1/2.2)).^2.5;
d.gamma(:,2) = (d.gamma(:,2).^(1/2.2)).^2.7;
d.gamma(:,3) = (d.gamma(:,3).^(1/2.2)).^2.3;

% COLOR MATCHING FUNCTIONS
S = [380 4 101]; % weird convention used by Brainard lab for defining wavelengths
load T_xyz1931; % load color matching functions
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale
wave = S(1):S(2):S(1)+S(2)*(S(3)-1); % define wavelength vector
% DEFOCUSES TO LOOK AT
Dall = -humanWaveDefocus(wave);
% WAVELENGTHS TO LOOK AT
% wvAll = humanWaveDefocusInvert(-Dall);

% PARAMETERS OF WAVEFRONT ANALYSIS
PARAMS.PixelDimension = 512;% size of pupil aperture field in pixels (this defines the resolution of the calculation)
PARAMS.PupilSize = 7; %default values - will be replaced depending on choices below
PARAMS.PupilFieldSize =42; %default values - will be replaced depending on choices below
PARAMS.PupilFitSize = 7; %default values - will be replaced depending on choices below
PARAMS.ImagingWavelength = 0.55;% imaging wavelength in microns
PARAMS.WavefrontResolution = 53;% increase to enhance the display of the wavefront (doesn't affect calculation)

%%

if subjNum==3
    subjName = 'S13-OD';
    blockNums = 12:17;
    trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
    % blockNums = [2 3];
    % trialNums = [[1:20]' [1:20]']; 
    nTrialTotal = 216;
elseif subjNum==10
    subjName = 'S20-OD';
    blockNums = 3:8;
    trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
    % blockNums = [2 3];
    % trialNums = [[1:20]' [1:20]'];     
    nTrialTotal = 216;
end

for l = 3 % LOOP OVER BLOCK
    for k = 1:36 % LOOP OVER TRIAL
        % LOADING DATA
        blockNumInd = l;
        blockNumTmp = blockNums(blockNumInd);
        trialNumTmp = trialNums(k,blockNumInd);
        
        AFCp = ARCloadFileBVAMS(subjNumEncode,blockNumTmp); % LOAD BVAMS DATA
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumTmp,0);

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS = struct;
        PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(:,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
        indBad = c(:,4)==0;
        meanC = mean(c(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS
        
        % STORE COLORS FOR FIRST AND SECOND STIMULI
        rgb00 = AFCp.rgb100(trialNumTmp,:);

        % recreate stimulus
        rVal = rgb00(1,1);
        gVal = rgb00(1,2);
        bVal = rgb00(1,3);
        im = imread([stimPath '/word_image_01.png']);
        im = double(im);
        imPattern = squeeze(im(:,:,3));
        imPattern = [zeros([60 size(imPattern,2)]); imPattern; zeros([60 size(imPattern,2)])];
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
        
        %% Computing peak correlation for different wavelengths in focus
        
        Dall2 = -humanWaveDefocus(wave(1:101));

        parfor i = 1:length(Dall2)
            zCoeffs = [0 meanC(1:end-1)];
            wvfP = wvfCreate('calc wavelengths', wave, ...
                'measured wavelength', humanWaveDefocusInvert(-Dall2(i)), ...
                'zcoeffs', zCoeffs, 'measured pupil', PARAMS.PupilSize, ...
                'name', sprintf('human-%d', PARAMS.PupilSize),'spatial samples',size(im,2));
            wvfP.calcpupilMM = PARAMS.PupilSize;
            defocusFromLCA = max(abs([humanWaveDefocusS10(humanWaveDefocusInvert(-Dall2(i)),min(wave)) ...
                              humanWaveDefocusS10(humanWaveDefocusInvert(-Dall2(i)),max(wave))]));
            if defocusFromLCA<1
                wvfP.refSizeOfFieldMM = 12;
            else
                wvfP.refSizeOfFieldMM = 6;
            end
            wvfP = wvfSet(wvfP, 'zcoeff', 0, 'defocus');
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS10);
            
            % Convert to siData format as well as wavefront object
            [siPSFData, wvfP] = wvf2SiPsfARC(wvfP,'showBar',false,'nPSFSamples',size(im,2),'umPerSample',1.1512); 
            oi = wvf2oi(wvfP); % CONVERT TO OPTICS OBJECT
            oi.optics.OTF.OTF = siPSFData.otf;
            % oi = oiCreateARC('human',wave,Dall(i)); % create optics
            oi = oiCompute(oi, s); % compute optical image of stimulus
        
            % Create the coneMosaic object
            cMosaic = coneMosaic;

            % Set size to show relevant portion of scene
            cMosaic.setSizeToFOV(1 * sceneGet(s, 'fov'));

            % key line for computing absorptions
            absorptions = cMosaic.computeSingleFrame(oi, 'fullLMS', true);            
            
            % absorptions = absorptions(55:128,6:177,:);

            display(['Peak correlation loop ' num2str(i) ' block ' num2str(blockNumTmp) ' trial ' num2str(k)]);
            absorptions = single(absorptions);
            S = struct;
            S.absorptions = absorptions;
            fnameCone = ['subj' num2str(subjNum) 'block' num2str(blockNumTmp) 'stimulus' num2str(k) 'focusInd' num2str(i)];
            save([savePath num2str(subjNum) '/' fnameCone '.mat'],"-fromstruct",S);
        end
    end
end

% save('/Users/benjaminchin/Documents/ARchromaScraps/ARCmodelOutput2_5.mat','meanv00all','wvInFocus1all','rgb1all','defocusBasic');
