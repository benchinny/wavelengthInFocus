function ARCpsfGeneration(subjNum,dataPath)

% subjNum values for participants who passed screening: 1, 3, 5, 10, 16,
% 17, 18, 20

%% Initialize and clear
ieInit;

% WHETHER OR NOT TO PLOT STIMULUS
bPlotStim = false;
% WHETHER OR NOT TO SAVE POINT-SPREAD FUNCTIONS
bSave = true;

%% Set up display struct and build Ben's stimulus

% PATH TO SAVE
savePath = fullfile(dataPath,'data','psfs');

% DEFINE LIGHT WAVELENGTHS TO SAMPLE STIMULUS AT
wave = 380:4:780;

% SET UP DISPLAY PARAMETERS (COMMON TO ALL RETINAL IMAGE MODELNG FOR THIS 
% PROJECT)
d = ARCmodelDispSetup(dataPath,0);

% LOAD COLOR CONDITIONS, PUPIL SIZE, AND MEAN ABERRATIONS
[~,rgbAll,~,PupilSize,meanC,~,~,~] = ARCnlzLoadDefocusAbb(subjNum,dataPath);

rgb00 = unique(rgbAll,'rows'); % GET UNIQUE COLOR CONDITIONS

% STIM COLOR DOESN'T REALLY MATTER HERE--WE JUST WANT THE STIMULUS TO
% FIGURE OUT HOW TO SIZE THE PSF
k = 12;
% MAKE STIMULUS USING COMMON HELPER FUNCTION
[s, ~] = ARCmodelStimSetup(dataPath,subjNum,'accommodation',d,rgb00(k,:),bPlotStim);
% PEAK CORRELATION WILL BE COMPUTED FOR A RANGE OF 'WAVELENGTHS IN
% FOCUS'. THIS IS A DISTINCT VARIABLE FROM 'wave', WHICH IS THE
% WAVELENGTHS OVER WHICH THE STIMULUS IS DEFINED. I HAVE MADE THEM THE
% SAME, BUT THEY DON'T NECESSARILY HAVE TO BE (E.G. IF YOU JUST WANT TO
% LOOK AT WHAT HAPPENS WHEN A SPECIFIC WAVELENGTH IS IN FOCUS)
wave2 = 380:4:780;

parfor i = 1:length(wave2) % LOOP OVER WAVELENGTHS IN FOCUS
    % % IF YOU WANT TO COMPUTE FOR DIFFRACTION-LIMITED SCENARIO 
    % zCoeffs = [0 zeros(size(meanC(1:end-1)))];

    % REFORMAT COEFFICIENTS (WAVEFRONT SENSOR LEAVES OUT THE 'PISTON'
    % ZERNIKE TERM, WHICH IS ALWAYS 0)
    zCoeffs = [0 meanC(1:end-1)];
    % NUMBER OF SPATIAL SAMPLES IN X AND Y (NOTE: NOT ROWS AND COLUMNS!)
    spatialSamplesXY = [size(s.data.photons,2) size(s.data.photons,1)];
    % FIX DEFOCUS TERM TO 0 BECAUSE WE ARE MODELING WAVELENGTH IN FOCUS
    defocusSet = 0;
    % CALL HELPER FUNCTION FOR SETTING UP OPTICS STRUCT SPECIFIC TO
    % MODELING FOR THIS PROJECT
    [oi, siPSFData] = ARCmodelOpticsSetup(subjNum,zCoeffs,wave,wave2(i),PupilSize,spatialSamplesXY,defocusSet);
    
    display(['Peak correlation loop ' num2str(i)]);

    % PLACE PSFS AND OTFS IN STRUCT
    S = struct;
    % STORE PSF
    S.psf = single(siPSFData.psf(134:184,134:184,:));
    % INITIALIZE MATRIX FOR STORING OTFS
    OTFzeroCentered = [];
    for j = 1:size(siPSFData.psf,3) % LOOP OVER WAVELENGTHS
        % STORE-ZERO-CENTERED OTF
        OTFzeroCentered(:,:,j) = ifftshift(squeeze(oi.optics.OTF.OTF(:,:,j)));
    end
    % JUST GRAB THE FIRST 60 FREQUENCIES
    S.otf = single(OTFzeroCentered(116:175,130:189,:));    
    % NAME FOR SAVING
    fnameCone = ['subj' num2str(subjNum) 'PSFfocusInd' num2str(i)];
    % UNCOMMENT LINE BELOW TO SAVE NEW CONE IMAGES
    if bSave
       save(fullfile(savePath,['S' num2str(subjNum)],[fnameCone '.mat']),"-fromstruct",S);
    end
end

end
