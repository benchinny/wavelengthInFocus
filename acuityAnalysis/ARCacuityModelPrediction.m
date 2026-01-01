function [dprimeMetric, dprime, dprimeCI] = ARCacuityModelPrediction(subjNum,LMOrChrom,dataPath)

% function for generating model predictions for acuity data.
%
% subjNum: subject number
% LMOrChrom: specify which model we are using to generate predictions.
%             'LM': luminance model.
%             'Chrom': best chromatic model
% dataPath  : path to folder where data lives
%
% dprimeMetric: d-prime predictions
% dprime      : empirical d-primes
% dprimeCI    : confidence intervals on empirical d-primes

%% Initialize and clear
ieInit;

% SET UP FOLDERS FOR SAVING AND LOADING MODEL PREDICTIONS FROM EXP1
saveFolder = fullfile(dataPath,'data','acuityModeling');
folderExp1 = fullfile(dataPath,'data','PresavedFigureData');
helperFolder = fullfile(dataPath,'data','helperFiles');
bPLOT = false;

if strcmp(LMOrChrom,'Chrom') % IF GENERATING CHROMATIC PREDICTIONS
    % LOAD BLUE-YELLOW PREDICTIONS
    load(fullfile(folderExp1,'wvMeanAndPredLMS.mat'),'dfPredPurpleAll','aicAll');
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
elseif strcmp(LMOrChrom,'LM') % IF USING LUMINANCE MODEL PREDICTIONS
    % LOAD THE LUMINANCE MODEL PREDICTIONS
    load(fullfile(folderExp1,'wvMeanAndPredLM.mat'),'dfPredPurpleAll','aicAll');
    modelPrediction875nmPurpleAt2pt5all = dfPredPurpleAll;
    % FOR SAVING FILE
    savePredName = 'acuityModelingPredictionLMS';
else
    error('LMOrChrom variable needs to either be the string ''LM'' or ''Chrom'' ');
end

% SET UP DISPLAY PARAMETERS (COMMON TO ALL RETINAL IMAGE MODELNG FOR THIS 
% PROJECT)
d = ARCmodelDispSetup(dataPath,0);

% SET UP STRUCTS FOR STIMULI
[s1, s2] = ARCmodelStimSetup(dataPath,subjNum,'acuity',d,[0.56 0 1],bPLOT);

% COLOR MATCHING FUNCTIONS
S = [380 4 101]; % weird convention used by Brainard lab for specifying wavelengths
% LOAD COLOR MATCHING FUNCTIONS--SECOND ROW IS V(LAMBDA)!
load T_xyz1931; 
% LOAD STANDARD LENS TRANSMITTANCE FUNCTION. WE WILL NEED TO FACTOR THIS OUT OF
% V-LAMBDA TO AVOID APPLYING TRANSMITTANCE TWICE.
load(fullfile(helperFolder,'transmittance.mat'));
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale
wave = S(1):S(2):S(1)+S(2)*(S(3)-1); % define wavelength vector
% FACTOR OUT TRANSMITTANCE FROM STANDARD V-LAMBDA
VlambdaNoTransmittance = T_sensorXYZ(2,:)./squeeze(transmittance)';

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
    % NUMBER OF SPATIAL SAMPLES IN X AND Y (NOTE: NOT ROWS AND COLUMNS!)
    spatialSamplesXY = [size(s1.data.photons,2) size(s1.data.photons,1)];    
    % SETTING DEFOCUS ACCORDING TO ACUITY STIMULUS DISTANCE
    defocusSet = -defocusForStim(i)*defocusScaleFactor;
    % SET UP OPTICS STRUCT
    oi = ARCmodelOpticsSetup(subjNum,zCoeffs,wave,875,PupilSize,spatialSamplesXY,defocusSet);
    
    oig1 = oiCompute(oi, s1); % compute optical image of stimulus
    oig2 = oiCompute(oi, s2); % compute optical image of stimulus
    
    % CONVERT FROM PHOTONS TO ENERGY TO LUMINANCE. FIRST, FORMAT TO XW
    % IMAGE FOR ISETBIO EXPECTATIONS
    photonsXW1 = RGB2XWFormat(oig1.data.photons);
    % CONVERT FROM QUANTA TO ENERGY
    energyXW1 = Quanta2Energy(wave,photonsXW1);
    % APPLY V LAMBDA AND RESHAPE TO XY IMAGE
    lumImgXW1 = sum(bsxfun(@times,energyXW1,VlambdaNoTransmittance),2);
    lumImgXY1 = reshape(lumImgXW1,[size(oig1.data.photons,1) size(oig1.data.photons,2)]);

    % CONVERT FROM PHOTONS TO ENERGY TO LUMINANCE. FIRST, FORMAT TO XW
    % IMAGE FOR ISETBIO EXPECTATIONS
    photonsXW2 = RGB2XWFormat(oig2.data.photons);
    % CONVERT FROM QUANTA TO ENERGY
    energyXW2 = Quanta2Energy(wave,photonsXW2);
    % APPLY V LAMBDA
    lumImgXW2 = sum(bsxfun(@times,energyXW2,VlambdaNoTransmittance),2);
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
