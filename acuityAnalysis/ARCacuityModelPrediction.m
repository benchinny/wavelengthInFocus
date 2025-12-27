function [dprimeMetric, dprime, dprimeCI] = ARCacuityModelPrediction(subjNum,LumOrChrom,dataPath)

% function for generating model predictions for acuity data.
%
% subjNum: subject number
% LumOrChrom: specify which model we are using to generate predictions.
%             'Lum': luminance model.
%             'Chrom': best chromatic model
% dataPath  : path to folder where data lives
%
% dprimeMetric: d-prime predictions
% dprime      : empirical d-primes
% dprimeCI    : confidence intervals on empirical d-primes

% FOR THIS FUNCTION TO GIVE CORRECT OUTPUT, MAKE SURE LENS TRANSMITTANCE 
% IN ISETBIO IS SET TO 1 EVERYWHERE! TO DO THIS, GO INTO THE FUNCTION 
% 'oiCalculateIrradiance' AND ADD THE LINE
% 'transmittance = ones(size(transmittance));' RIGHT BEFORE LINE 87. 

%% Initialize and clear
ieInit;

% SET UP FOLDERS FOR SAVING AND LOADING MODEL PREDICTIONS FROM EXP1
saveFolder = fullfile(dataPath,'data','acuityModeling');
folderExp1 = fullfile(dataPath,'data','PresavedFigureData');
bPLOT = false;

if strcmp(LumOrChrom,'Chrom') % IF GENERATING CHROMATIC PREDICTIONS
    % LOAD BLUE-YELLOW PREDICTIONS
    load(fullfile(folderExp1,'wvMeanAndPredDonutx2.mat'),'dfPredPurpleAll','aicAll');
    dfPredPurpleBYAll = dfPredPurpleAll;
    aicBYall = aicAll;
    % LOAD RED-GREEN PREDICTIONS
    load(fullfile(folderExp1,'wvMeanAndPredLminusM.mat'),'dfPredPurpleAll','aicAll');
    dfPredPurpleRGAll = dfPredPurpleAll;
    aicRGall = aicAll;    
    % USE AIC TO DETERMINE BETTER CHROMATIC MODEL FOR EACH SUBJECT
    modelPrediction875nmPurpleAt2pt5all = zeros(size(dfPredPurpleBYAll));
    indBYbetter = aicBYall<aicRGall;
    % STORE BEST CHROMATIC MODEL PREDICTIONS
    modelPrediction875nmPurpleAt2pt5all(indBYbetter) = dfPredPurpleBYAll(indBYbetter);
    modelPrediction875nmPurpleAt2pt5all(~indBYbetter) = dfPredPurpleRGAll(~indBYbetter);
    % FOR SAVING FILE
    savePredName = 'acuityModelingPredictionS';
elseif strcmp(LumOrChrom,'Lum') % IF USING LUMINANCE MODEL PREDICTIONS
    % LOAD THE LUMINANCE MODEL PREDICTIONS
    load(fullfile(folderExp1,'wvMeanAndPredLM.mat'),'dfPredPurpleAll','aicAll');
    modelPrediction875nmPurpleAt2pt5all = dfPredPurpleAll;
    % FOR SAVING FILE
    savePredName = 'acuityModelingPredictionLumS';
else
    error('LumOrChrom variable needs to either be the string ''Lum'' or ''Chrom'' ');
end

% CALIBRATION DATA FILE PATH
calPath = fullfile(dataPath,'BVAMS_calibration_files','Ben_calibration_July_6_2024');

% SET UP DISPLAY PARAMETERS (COMMON TO ALL RETINAL IMAGE MODELNG FOR THIS 
% PROJECT)
d = ARCmodelDispSetup(calPath);

% RGB VALUES OF ACUITY STIMULUS IN ACUITY TASK
rVal = 0.56;
bVal = 1.00;
gVal = 0.00;

% GET THRESHOLDS FROM CONTRAST CALIBRATION EXPERIMENT
thresholds = ARCnlz_contrastThresholds(subjNum,0,dataPath);

% GABOR PARAMETERS
frqCpd = 15;
contrast = thresholds(4);
rgbAll = [rVal gVal bVal];
% GAMMA (BASED ON BVAMS CALIBRATION)
gammaR = 2.5;
gammaG = 2.7;
gammaB = 2.3;

% STIMULUS PARAMETERS 
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
               [rgbAll(1,1)^gammaR rgbAll(1,2)^gammaG rgbAll(1,3)^gammaB],1,1,0,0);

% MAKE SURE TO GAMMA CORRECT STIMULUS FOR MODELING
acuStimOrig1(:,:,1) = acuStimOrig1(:,:,1).^(1/gammaR);
acuStimOrig1(:,:,2) = acuStimOrig1(:,:,2).^(1/gammaG);
acuStimOrig1(:,:,3) = acuStimOrig1(:,:,3).^(1/gammaB);
I1 = acuStimOrig1;

% MODEL ACUITY STIMULUS FOR -15 DEG ORIENTATION
acuStimOrig2 = ARC2Dgabor(stimPositionsX,[],x0,y0,frqCpdAll, ...
               contrastAll,orientation2,phs,sigmaX,sigmaY, ...
               [rgbAll(1,1)^gammaR rgbAll(1,2)^gammaG rgbAll(1,3)^gammaB],1,1,0,0);

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

if bPLOT
    figure; 
    set(gcf,'Position',[289 428 1056 420]);
    subplot(1,3,1);
    plot(d.wave,d.spd(:,1),'r','LineWidth',1.5); hold on;
    plot(d.wave,d.spd(:,2),'g','LineWidth',1.5);
    plot(d.wave,d.spd(:,3),'b','LineWidth',1.5);
    axis square;
    formatFigure('Wavelength (\lambda)','Radiance');
    subplot(1,3,2);
    imagesc(I1);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    axis square;
    set(gca,'FontSize',15);
    title('Original');
    subplot(1,3,3);
    plot(s1.spectrum.wave,squeeze(s1.data.photons(130,130,:)),'-k','LineWidth',1);
    formatFigure('Wavelength (\lambda)','Photons');
    axis square;
end

% COLOR MATCHING FUNCTIONS
S = [380 4 101]; % weird convention used by Brainard lab for specifying wavelengths
load T_xyz1931; % load color matching functions
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale
wave = S(1):S(2):S(1)+S(2)*(S(3)-1); % define wavelength vector

% GET ZERNIKE COEFFICIENTS FOR PARTICIPANT
[~, ~, ~, PupilSize, meanC] = ARCnlzLoadDefocusAbb(subjNum,dataPath);

dprimeMetric = []; % INITIALIZING VECTOR FOR STORING D-PRIME VALUES
defocusScaleFactor = 0.5774; % FOR 4MM PUPIL SIZE

subjNumAll = [1 3 5 10 16 17 18 20]; % LIST OF ALL SUBJECT IDS
% GRAB PREDICTION OF MODEL AT 2.5D FOR THAT PARTICULAR SUBJECT (BASED ON
% PRE-LOADED DATA EARLY IN THIS FUNCTION)
modelPrediction875nmPurpleAt2pt5 = -modelPrediction875nmPurpleAt2pt5all(subjNum==subjNumAll);

% MODEL ACUITY STIMULUS AT DIFFERENT DISTANCES
defocusForStim = [0.6:0.1:4.4]-modelPrediction875nmPurpleAt2pt5;

% LOAD PRE-SAVED LCA PARAMETERS
load(fullfile(dataPath,'data','PresavedFigureData','LCAparams.mat'),'q1bestAll','q2bestAll','q3bestAll');
q1 = q1bestAll(subjNum==subjNumAll);
q2 = q2bestAll(subjNum==subjNumAll);
q3 = q3bestAll(subjNum==subjNumAll);

% CONVERT TO WAVELENGTH-IN-FOCUS
wvInFocusForStim = humanWaveDefocusInvertParameterizedARC(875,-defocusForStim,q1,q2,q3);

parfor i = 1:length(defocusForStim)
    % FORMAT ACCORDING TO WHAT ISETBIO EXPECTS (WAVEFRONT SENSOR LEAVES OUT
    % PISTON TERM, ISETBIO DOESN'T)
    zCoeffs = [0 meanC(1:end-1)];
    wvfP = wvfCreate('calc wavelengths', wave, ...
        'measured wavelength', 875, ...
        'zcoeffs', zCoeffs, 'measured pupil', PupilSize, ...
        'name', sprintf('human-%d', PupilSize),'spatial samples',size(I1,2));
    % MAKE SURE THIS VARIABLE IS SET TO THE ACTUAL PUPIL SIZE
    wvfP.calcpupilMM = PupilSize;
    % THE AREA OF THE PUPIL THE PSF IS CALCULATED FROM
    wvfP.refSizeOfFieldMM = 6;    
    % SET ZERNIKE COEFFICIENT ACCORDING TO STIMULUS DISTANCE
    wvfP = wvfSet(wvfP, 'zcoeff', -defocusForStim(i)*defocusScaleFactor, 'defocus');
    % SET CUSTOM LCA FUNCTION PER SUBJECT--ISETBIO WANTS IT TO BE SET
    % IN A PARTICULAR FORMAT
    if subjNum==10
        wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS10);
    elseif subjNum==3    
        wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS3);
    elseif subjNum==5     
        wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS5);
    elseif subjNum==1      
        wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS1);
    elseif subjNum==16      
        wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS16);
    elseif subjNum==17   
        wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS17); 
    elseif subjNum==18      
        wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS18); 
    elseif subjNum==20
        wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS20);         
    end
    
    % MAKE POINT-SPREAD FUNCTION (siPSFData) AND WAVEFRONT STRUCT
    [siPSFData, wvfP] = wvf2SiPsfARC(wvfP,'showBar',false,'nPSFSamples',size(I1,2),'umPerSample',1.1512); % 1.1512
    oi = wvf2oi(wvfP); % CONVERT TO OPTICS OBJECT
    % NEED TO REMOVE PADDED ZEROS FROM PSF TO MAKE SAME SIZE AS
    % STIMULUS IMAGE. THE LINES BELOW IDENTIFY THE 'GOOD INDICES', I.E.
    % THE INDICES THAT AREN'T THE PADDED ZEROS TO BE REMOVED
    paddingXCpsf = round((size(siPSFData.psf,2)-size(s1.data.photons,2))/2);
    paddingYRpsf = round((size(siPSFData.psf,1)-size(s1.data.photons,1))/2); 
    indNotPadded = {(paddingYRpsf+1):(size(siPSFData.psf,1)-paddingYRpsf) ...
                    (paddingXCpsf+1):(size(siPSFData.psf,2)-paddingXCpsf)}; 
    % REPLACE OTF FIELD WITH A NEW OTF CALCULATED FROM THE PSFS WE
    % JUST CALCULATED. COULDN'T FIGURE OUT HOW TO MAKE ISETBIO DO
    % THIS AUTOMATICALLY, SO I DID IT MANUALLY (AND CHECKED THE
    % OUTPUT!) 
    oi.optics.OTF = []; % INITIALIZE ARRAY FOR STORING OTFS
    for j = 1:size(siPSFData.psf,3) % LOOP OVER WAVELENGTHS
        % NOTE THAT WE APPLY fftshift TO THE PSF SO THAT ITS CENTER IS
        % IN INDEX (1,1) OF THE IMAGE (TOP LEFT). fft2 EXPECTS THE
        % SIGNAL ORIGIN TO BE IN THIS LOCATION.        
        oi.optics.OTF.OTF(:,:,j) = fft2(fftshift(squeeze(siPSFData.psf(indNotPadded{1},indNotPadded{2},j))));
    end
    % THESE TWO STEPS WILL YIELD INCORRECT OUTPUT IF LENS TRANSMITTANCE IS
    % NOT SET TO 1 THROUGHOUT (SEE COMMENTS AT TOP OF FUNCTION)
    oig1 = oiCompute(oi, s1); % compute optical image of stimulus
    oig2 = oiCompute(oi, s2); % compute optical image of stimulus
    
    % CONVERT FROM PHOTONS TO ENERGY TO LUMINANCE. FIRST, FORMAT TO XW
    % IMAGE FOR ISETBIO EXPECTATIONS
    photonsXW1 = RGB2XWFormat(oig1.data.photons);
    % CONVERT FROM QUANTA TO ENERGY
    energyXW1 = Quanta2Energy(wave,photonsXW1);
    % APPLY V LAMBDA AND RESHAPE TO XY IMAGE
    lumImgXW1 = sum(bsxfun(@times,energyXW1,squeeze(T_sensorXYZ(2,:))),2);
    lumImgXY1 = reshape(lumImgXW1,[size(oig1.data.photons,1) size(oig1.data.photons,2)]);

    % CONVERT FROM PHOTONS TO ENERGY TO LUMINANCE. FIRST, FORMAT TO XW
    % IMAGE FOR ISETBIO EXPECTATIONS
    photonsXW2 = RGB2XWFormat(oig2.data.photons);
    % CONVERT FROM QUANTA TO ENERGY
    energyXW2 = Quanta2Energy(wave,photonsXW2);
    % APPLY V LAMBDA
    lumImgXW2 = sum(bsxfun(@times,energyXW2,squeeze(T_sensorXYZ(2,:))),2);
    % NO NEED TO PLOT STIMULUS IN BOTH ORIENTATIONS UNLESS YOU WANT TO
    % lumImgXY2 = reshape(lumImgXW2,[size(oig2.data.photons,1) size(oig2.data.photons,2)]);
    
    % CONVERT TO D-PRIME ACCORDING TO GEISLER (1989) EQUATION
    dprimeMetricDenom = sqrt(sum(sum((lumImgXW2+lumImgXW1).*log(lumImgXW2./lumImgXW1).^2)));
    dprimeMetricNumer = sum(sum((lumImgXW2-lumImgXW1).*log(lumImgXW2./lumImgXW1)));
    dprimeMetric(i) = dprimeMetricNumer./dprimeMetricDenom;
    
    if bPLOT
        figure; 
        set(gcf,'Position',[359 443 975 420]);
        subplot(1,2,1);
        imagesc(lumImgXY1); 
        axis square; 
        colormap gray;
        title(['Defocus = ' num2str(defocusForStim(i)) ', x-correlation = ' num2str(dprimeMetric(i))]);
        subplot(1,2,2);
        imagesc(fftshift(ifft2(squeeze(oig1.optics.OTF.OTF(:,:,61))))); 
        axis square; 
        colormap gray;
    end
    display(['D-prime iteration ' num2str(i)]);
end

% EMPIRICAL PERFORMANCE INCLUDING CALCULATED D-PRIME VALUES
[unqFocDst,PC,PCci,dprime,dprimeCI,PCfit,dprimeFitAll,PCfitSupport] = ARCacuityAnalyzeDataOnly(subjNum,0,dataPath);

if bPLOT
    figure;
    set(gcf,'Position',[342 460 1052 440]);
    subplot(1,2,1);
    hold on;
    plot(defocusForStim+modelPrediction875nmPurpleAt2pt5,dprimeMetric,'-','Color',[0.56 0 1],'LineWidth',1);
    scaleFac = 0.816;
    dprimeScale = max(dprime(:)./max(dprimeMetric));
    errorbar(2.5+unqFocDst.*scaleFac,dprime./dprimeScale,(dprime-dprimeCI(1,:))./dprimeScale,(dprimeCI(2,:)-dprime)./dprimeScale,'o','Color',[0.56 0 1],'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
    xlabel('Distance');
    ylabel('D-prime metric');
    set(gca,'FontSize',15);
    title(['Mean defocus at 875nm = ' num2str(modelPrediction875nmPurpleAt2pt5,3) 'D']);
    subplot(1,2,2);
    plot(wvInFocusForStim,dprimeMetric,'k-','LineWidth',1);
    xlabel('Wavelength in focus (nm)');
    ylabel('D-prime metric');
    set(gca,'FontSize',15);
end

save(fullfile(saveFolder,[savePredName num2str(subjNum)]),'dprimeMetric','defocusForStim', ...
    'modelPrediction875nmPurpleAt2pt5','dprime','dprimeCI','unqFocDst','wvInFocusForStim');

end
