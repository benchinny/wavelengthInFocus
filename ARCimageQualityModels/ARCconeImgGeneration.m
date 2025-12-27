function oi = ARCconeImgGeneration(subjNum,dataPath)

% subjNum values for participants who passed screening: 1, 3, 5, 10, 16,
% 17, 18, 20

% MAKE SURE LENS TRANSMITTANCE IN ISETBIO IS DEFAULT! TO DO
% THIS, GO INTO THE FUNCTION 'oiCalculateIrradiance' AND MAKE SURE THERE IS
% NO LINE THAT SAYS 'transmittance = ones(size(transmittance));' RIGHT 
% BEFORE LINE 87. 

%% Initialize and clear
ieInit;

% WHETHER OR NOT TO PLOT STIMULUS
bPlotStim = false;
bSave = true;

%% Set up display struct and build Ben's stimulus

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
        % CREATE ISETBIO WAVEFRONT OBJECT
        wvfP = wvfCreate('calc wavelengths', wave, ...
            'measured wavelength', wave2(i), ...
            'zcoeffs', zCoeffs, 'measured pupil', PupilSize, ...
            'name', sprintf('human-%d', PupilSize),'spatial samples',size(s.data.photons,2));
        % MAKE SURE THE calcpupilMM FIELD MATCHES THE ACTUAL PUPIL SIZE IN
        % THE EXPERIMENT
        wvfP.calcpupilMM = PupilSize;
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
        [siPSFData, wvfP] = wvf2SiPsfARC(wvfP,'showBar',false,'nPSFSamples',size(s.data.photons,2),'umPerSample',1.1512); 
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
            % NOTE THAT WE APPLY fftshift TO THE PSF SO THAT ITS CENTER IS
            % IN INDEX (1,1) OF THE IMAGE (TOP LEFT). fft2 EXPECTS THE
            % SIGNAL ORIGIN TO BE IN THIS LOCATION.
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
        if bSave
           save(fullfile(savePath,['S' num2str(subjNum)],[fnameCone '.mat']),"-fromstruct",S);
        end
    end
end

end
