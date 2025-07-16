function oi = ARCimgQualityAnalysisConesMeanZfinch(subjNum,primaryCell)

% NOTE SUBJECT NUMBER CONVENTION: SUBTRACT 10 FROM subjNum TO GET ACTUAL
% SUBJECT NUMBER. subjNum VALUES <=10 WERE INTENTIONALLY NOT USED FOR
% ACTUAL PARTICIPANTS. NOTE ALSO THAT PARTICIPANTS WHO DID NOT PASS
% SCREENING OR HAD TO BE EXCLUDED FROM THE ACTUAL ANALYSIS ARE STILL
% INCLUDED IN THIS FUNCTION. 

% subjNum values for participants who passed screening: 11, 13, 15, 20, 26,
% 27, 28, 30. That is, subjects S1, S3, S5, S10, S16, S17, S18, S20. 

%% Initialize and clear
ieInit;

%% Set up display struct and build Ben's stimulus

subjNumEncode = subjNum+10;
bPlotCoefficients = false;

% Setting up display properties
d = displayCreate('OLED-Samsung');
d = displaySet(d, 'name', 'my display');
d = displaySet(d,'ViewingDistance',1); % simulated screen distance
d = displaySet(d,'dpi',378); % simulated screen distance

bUseBVAMScal = 1; % if using BVAMS calibration data

if strcmp(getenv('USER'),'benjaminchin')
    calPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/';
    stimPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/';
    savePath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImagesFinch/S';
end

if strcmp(getenv('USER'),'benchin')
    calPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/';
    stimPath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/Meeting_April_23/';
    savePath = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneImagesFinch/S';
end

% COLOR MATCHING FUNCTIONS
S = [380 4 101]; % weird convention used by Brainard lab for defining wavelengths
load T_xyz1931; % load color matching functions
T_sensorXYZ = 683*SplineCmf(S_xyz1931,T_xyz1931,S); % interpolate and scale
wave = S(1):S(2):S(1)+S(2)*(S(3)-1); % define wavelength vector
% DEFOCUSES TO LOOK AT
Dall = -humanWaveDefocus(wave);

% primaryCell = {'red' 'green'};
if bUseBVAMScal
    load([calPath 'Finch_et_al_primaries_' primaryCell{1} '.mat']);
    spd1 = zeros(size(wave));
    indSpd1 = wave>min(waveInterp) & wave<max(waveInterp);
    spd1(indSpd1) = interp1(waveInterp,powerInterp,wave(indSpd1));
    load([calPath 'Finch_et_al_primaries_' primaryCell{2} '.mat']);
    spd2 = zeros(size(wave));
    indSpd2 = wave>min(waveInterp) & wave<max(waveInterp);
    spd2(indSpd2) = interp1(waveInterp,powerInterp,wave(indSpd2));    
    d.spd = zeros([length(spd1) 3]);
    d.spd(:,1) = spd1;
    d.spd(:,3) = spd2;
end
d.gamma(:,1) = (d.gamma(:,1).^(1/2.2)).^2.2;
d.gamma(:,2) = (d.gamma(:,2).^(1/2.2)).^2.2;
d.gamma(:,3) = (d.gamma(:,3).^(1/2.2)).^2.2;

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
elseif subjNum==1
   blockNums = 11:16;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
   subjName = ['S' num2str(subjNum+10) '-OD'];
   nTrialTotal = 216;
elseif subjNum==5
   blockNums = 3:8;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
   subjName = ['S' num2str(subjNum+10) '-OD'];
   nTrialTotal = 216;   
elseif subjNum==9
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]']; 
   subjName = ['S' num2str(subjNum+10) '-OD'];
   nTrialTotal = 216;      
elseif subjNum==16
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+10) '-OD'];
   nTrialTotal = 216;
elseif subjNum==17
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+10) '-OD'];
   nTrialTotal = 216;
elseif subjNum==18
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+10) '-OD'];
   nTrialTotal = 216; 
elseif subjNum==20
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+10) '-OD'];
   nTrialTotal = 216;
elseif subjNum==21
   blockNums = 2:7;
   trialNums = [[1:36]' [1:36]' [1:36]' [1:36]' [1:36]' [1:36]'];
   subjName = ['S' num2str(subjNum+9) '-OD'];
   nTrialTotal = 216;   
   subjNumEncode = subjNum+9;
end

%%

lumPropRB = [0.000 1.000; ...
             0.125 0.875; ...
             0.250 0.750; ...
             0.375 0.625; ...
             0.500 0.500; ...
             0.625 0.375; ...
             0.750 0.250; ...
             0.875 0.125; ...
             1.000 0.000];
lumB = (T_sensorXYZ(2,:)*d.spd(:,3));
lumR = (T_sensorXYZ(2,:)*d.spd(:,1));

if lumB>lumR
    lumPropRB(:,2) = lumPropRB(:,2).*(lumR./lumB);
else
    lumPropRB(:,1) = lumPropRB(:,1).*(lumB./lumR);
end
rgb00 = zeros([9 3]);
rgb00(:,1) = lumPropRB(:,1).^(1/2.2);
rgb00(:,3) = lumPropRB(:,2).^(1/2.2);

cAll = [];
optDistAll = [];
rgbAll = [];

for l = 1:length(blockNums) % LOOP OVER BLOCK
    blockNumInd = l;
    blockNumTmp = blockNums(blockNumInd);
    AFCp = ARCloadFileBVAMS(subjNumEncode,blockNumTmp); % LOAD BVAMS DATA
    rgbAll = [rgbAll; AFCp.rgb100];
    for k = 1:36 % LOOP OVER TRIAL
        trialNumTmp = trialNums(k,blockNumInd);
        
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumTmp,0);

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS = struct;
        PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
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
% rgb00 = unique(rgbAll,'rows');
if subjNum==21
    meanC = [0 0 -0.036294 0.8333 -0.16442 -0.0063637 0.0063415 -0.0086154 0.02253 0.0094638 0.004479 0.01336 -0.01776 0.01666 -0.0044095 -0.0024943 0.0044616 0.0058616 -0.0024607 0.0060126 -0.0040849 5.2909e-07 0.001748 -0.0061759 0.0028174 -0.0052323 -0.0065851 0.00098404 0.0031324 0.0013516 -0.0012985 -0.0046191 -2.1413e-05 0.0011138 0.00079841 0.0043567 0.00079581 -0.0010823 -0.0017864 0.0029658 -0.00098402 0.0022702 0.0020892 -0.0023476 -0.0010185 0.00010112 -0.00084987 -0.00049893 -9.0462e-05 0.0018192 -0.00023171 -0.00075573 4.7396e-05 -0.00025469 -0.0019559 0.00081874 -3.6993e-05 -0.00042477 0.00025839 -0.00048126 0.00033523 -0.0007869 -0.00015903 0.00078828 -0.0015525];
end

%% PLOTTING EACH COEFFICIENT VS DEFOCUS 

if bPlotCoefficients
    figure;
    set(gcf,'Position',[163 358 1304 588]);
    for i = 1:10
        subplot(2,5,i);
        hist(cAll(~indBad,i+2),linspace(-0.5,2,101));
        axis square;
    end
    
    figure;
    set(gcf,'Position',[117 149 1393 773]);
    subplot(2,1,1);
    hold on;
    for i = 1:65
        bar(i,mean(cAll(~indBad,i)));
    end
    formatFigure('Zernike polynomial #','Mean coefficient');
    subplot(2,1,2);
    hold on;
    for i = 1:65
        bar(i,std(cAll(~indBad,i)));
    end
    formatFigure('Zernike polynomial #','Std dev coefficient');
    
    coeffs2compare = [3 5:13];
    figure;
    set(gcf,'Position',[117 149 1393 773]);
    for i = 1:length(coeffs2compare)
       subplot(2,5,i);
       plot(cAll(~indBad,4),cAll(~indBad,coeffs2compare(i)),'k.');
       % ylim([-0.25 0.25]);
       axis square;
       if i==1
          formatFigure('Defocus coefficient',['Coefficient ' num2str(coeffs2compare(i))],'Astigmatism');
       elseif i==2
          formatFigure('Defocus coefficient',['Coefficient ' num2str(coeffs2compare(i))],'Astigmatism');
       elseif i==3
          formatFigure('Defocus coefficient',['Coefficient ' num2str(coeffs2compare(i))],'Trefoil 0');
       elseif i==4
          formatFigure('Defocus coefficient',['Coefficient ' num2str(coeffs2compare(i))],'Coma Y');  
       elseif i==5
          formatFigure('Defocus coefficient',['Coefficient ' num2str(coeffs2compare(i))],'Coma X');
       elseif i==6
          formatFigure('Defocus coefficient',['Coefficient ' num2str(coeffs2compare(i))],'Trefoil 30');
       elseif i==7 || i==8
          formatFigure('Defocus coefficient',['Coefficient ' num2str(coeffs2compare(i))],'');
       elseif i==9
          formatFigure('Defocus coefficient',['Coefficient ' num2str(coeffs2compare(i))],'Spherical aberration');      
       else
          formatFigure('Defocus coefficient',['Coefficient ' num2str(coeffs2compare(i))]);
       end
    end

    % PLOTTING EACH COEFFICIENT VS OPTICAL DISTANCE
    
    optDistUnq = unique(optDistAll);
    
    figure;
    set(gcf,'Position',[117 149 1393 773]);
    for i = 1:length(coeffs2compare)
       subplot(2,5,i);
       hold on;
       for j = 1:length(optDistUnq)
           indOptDist = abs(optDistAll-optDistUnq(j))<0.001;
           errorbar(optDistUnq(j),mean(cAll(~indBad & indOptDist,coeffs2compare(i))), ...
                    std(cAll(~indBad & indOptDist,coeffs2compare(i))),'ko', ...
                    'MarkerSize',10,'MarkerFaceColor','w','LineWidth',1);
       end
       axis square;
       xlim([1 4]);
       if i==1
          formatFigure('Optical distance (D)',['Coefficient ' num2str(coeffs2compare(i))],'Astigmatism');
       elseif i==2
          formatFigure('Optical distance (D)',['Coefficient ' num2str(coeffs2compare(i))],'Astigmatism');
       elseif i==3
          formatFigure('Optical distance (D)',['Coefficient ' num2str(coeffs2compare(i))],'Trefoil 0');
       elseif i==4
          formatFigure('Optical distance (D)',['Coefficient ' num2str(coeffs2compare(i))],'Coma Y');  
       elseif i==5
          formatFigure('Optical distance (D)',['Coefficient ' num2str(coeffs2compare(i))],'Coma X');
       elseif i==6
          formatFigure('Optical distance (D)',['Coefficient ' num2str(coeffs2compare(i))],'Trefoil 30');
       elseif i==7 || i==8
          formatFigure('Optical distance (D)',['Coefficient ' num2str(coeffs2compare(i))],'');
       elseif i==9
          formatFigure('Optical distance (D)',['Coefficient ' num2str(coeffs2compare(i))],'Spherical aberration');      
       else
          formatFigure('Optical distance (D)',['Coefficient ' num2str(coeffs2compare(i))]);
       end   
    end
end

%%

for k = 1:9 % LOOP OVER TRIAL
    % recreate stimulus
    rVal = rgb00(k,1);
    gVal = rgb00(k,2);
    bVal = rgb00(k,3);
    im = imread([stimPath '/Finch_et_al_stimulus.png']);
    im = imresize(im,[390 390]);
    im = double(im);
    imPattern = squeeze(im(:,:,3));
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
    
    wave2 = 380:4:780;

    parfor i = 1:length(wave2)
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
        elseif subjNum==21
            defocusFromLCA = max(abs([humanWaveDefocusS21(wave2(i),min(wave)) ...
                                      humanWaveDefocusS21(wave2(i),max(wave))]));  
            wvfP = wvfSet(wvfP, 'customlca', @humanWaveDefocusS21);    
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
        
        % absorptions = absorptions(55:128,6:177,:);

        display(['Peak correlation loop ' num2str(i) ' stimulus ' num2str(k)]);
        absorptions = single(absorptions);
        S = struct;
        S.absorptions = absorptions;
        fnameCone = ['subj' num2str(subjNum) 'stimulus' num2str(k) primaryCell{1} primaryCell{2} 'focusInd' num2str(i)];
        save([savePath num2str(subjNum) '/' fnameCone '.mat'],"-fromstruct",S);
    end
end

end
