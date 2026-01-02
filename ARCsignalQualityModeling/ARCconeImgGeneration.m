function ARCconeImgGeneration(subjNum,dataPath)

% function for generating cone images
%
% subjNum values for participants who passed screening: 1, 3, 5, 10, 16,
% 17, 18, 20

%% Initialize and clear
ieInit;

% WHETHER OR NOT TO PLOT STIMULUS
bPlotStim = false;
% WHETHER OR NOT TO SAVE CONE IMAGES
bSave = true;
% WHETHER OR NOT TO APPLY OPTICS
bApplyOptics = true;

%% Set up display struct and build stimulus

% PATH TO SAVE
savePath = fullfile(dataPath,'data','coneImages');

% DEFINE LIGHT WAVELENGTHS TO SAMPLE STIMULUS AT
wave = 380:4:780;

% SET UP DISPLAY PARAMETERS (COMMON TO ALL RETINAL IMAGE MODELNG FOR THIS 
% PROJECT)
d = ARCmodelDispSetup(dataPath,0);

% LOAD COLOR CONDITIONS, PUPIL SIZE, AND MEAN ABERRATIONS
[~,rgbAll,~,PupilSize,meanC,~,~,~] = ARCnlzLoadDefocusAbb(subjNum,dataPath);

rgb00 = unique(rgbAll,'rows'); % GET UNIQUE COLOR CONDITIONS

for k = 1:size(rgb00,1) % LOOP OVER COLOR CONDITIONS
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
        oi = ARCmodelOpticsSetup(subjNum,zCoeffs,wave,wave2(i),PupilSize,spatialSamplesXY,defocusSet);

        if bApplyOptics % IF APPLYING OPTICS
           % MAIN STEP: COMPUTE OPTICAL IMAGE OF STIMULUS
           oi = oiCompute(oi, s); 
           % NO EXTRA TAGS IN FILE NAME
           opticsNameTag = '';
        else % IF NOT APPLYING OPTICS
           % MAIN STEP: COMPUTE OPTICAL IMAGE OF STIMULUS
           oi = oiComputeNoOptics(oi,s);
           % TAG FILENAME TO INDICATE NO OPTICS WERE USED
           opticsNameTag = 'noOptics';
        end

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
        fnameCone = ['subj' num2str(subjNum) 'stimulus' num2str(k) 'focusInd' num2str(i) opticsNameTag];
        if bSave
           save(fullfile(savePath,['S' num2str(subjNum)],[fnameCone '.mat']),"-fromstruct",S);
        end
    end
end

end
