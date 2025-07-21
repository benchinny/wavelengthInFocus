%%220508 JnJ_AFC AFC9f include TCA correction. 
ex='RGB'; %ey=input('which eye? Right/Binc');
ey=input('which eye are we testing (type 1 for Right or 2 for Binocular)? ');
if ey==1; ey='Right'; elseif ey==2; ey='Binc'; end
filePath = 'H:\Shared drives\CIVO_BVAMS\data\ARC\';
vsIncrement = input(['Increment visit number? The current visit number is ' num2str(vs) ' (1 for yes, 0 for no)']);
vsIncrement = input(['Asking again: increment visit number? The current visit number is ' num2str(vs) ' (1 for yes, 0 for no)']);

if vsIncrement>=1
   vs = vs+1;
end

% expSubType = 'OLD';
% expSubType = 'STP';
% expSubType = 'SIN';
% expSubType = 'STEP';
expSubType = 'RGB';
blockType = 'random';

meanFocstmOptDstUnq = [1.5 2.0 2.5 3.0 3.5]';
meanFocstmOptDstUnq = meanFocstmOptDstUnq*1.25;

if strcmp(blockType,'fixed')
    nRepeats = 5;
    rgb1 = [0.555 0.245 0.533; ...
            0.555 0.245 0.533; ...
            0.312 0.418 0.533; ...
            0.312 0.418 0.533; ...
            0.312 0.245 1.00; ...   
            0.312 0.245 1.00];
    rgb2 = [0.555 0.245 0.533; ...
            0.555 0.245 0.533; ...
            0.312 0.418 0.533; ...
            0.312 0.418 0.533; ...
            0.312 0.245 1.00; ...            
            0.312 0.245 1.00];
           
    % TRIALS
    focStmOptDst = 1.25.*[1 -1 1 -1 1 -1]';
    meanFocstmOptDst = 2.5.*1.25.*ones([size(rgb1,1) 1]); 
    rgb1 = repmat(rgb1,[nRepeats 1]);
    rgb2 = repmat(rgb2,[nRepeats 1]);
    focStmOptDst = repmat(focStmOptDst,[nRepeats 1]);
    meanFocstmOptDst = repmat(meanFocstmOptDst,[nRepeats 1]);
    AFCv = randsample(1:length(focStmOptDst),length(focStmOptDst));
end

if strcmp(blockType,'random')
    % nCnd = 20; % NUMBER OF CONDITIONS
    % rVals = [0.56 0.42 0.32 0.24]'; % POSSIBLE RED PRIMARY VALUES
    % bVals = [1.00 0.73 0.53 0.39]'; % POSSIBLE BLUE PRIMARY VALUES
    % stepVals = 1.25.*[1.25 0.75 -0.75 -1.25]'; % POSSIBLE STEP VALUES
    % indCnd = reshape(randsample(1:length(rVals),nCnd*5,1),[nCnd 5]);
    % rgb1 = [rVals(indCnd(:,1)) zeros([nCnd 1]) bVals(indCnd(:,2))];
    % rgb2 = [rVals(indCnd(:,3)) zeros([nCnd 1]) bVals(indCnd(:,4))];
    % meanFocstmOptDst = meanFocstmOptDst.*ones([nCnd 1]);
    % focStmOptDst = stepVals(indCnd(:,5)).*ones([nCnd 1]);
    % AFCv = 1:length(focStmOptDst);
    nCnd = 20; % NUMBER OF CONDITIONS
    rVals = [0.555 0.492 0.416 0.312]'; % POSSIBLE RED PRIMARY VALUES
    gVals = [0.418 0.375 0.320 0.245]'; % POSSIBLE GREEN PRIMARY VALUES
    bVals = [1.00 0.877 0.730 0.533]'; % POSSIBLE BLUE PRIMARY VALUES
    % stepVals = 1.25.*[1.25 0.75 -0.75 -1.25]'; % POSSIBLE STEP VALUES
    % RANDOMLY PICKING COLORS
    indCnd = reshape(randsample(1:length(rVals),nCnd*7,1),[nCnd 7]);
    rgb1 = [rVals(indCnd(:,1)) gVals(indCnd(:,2)) bVals(indCnd(:,3))];
    rgb2 = [rVals(indCnd(:,4)) gVals(indCnd(:,5)) bVals(indCnd(:,6))];
    % RANDOMLY PICKING OPTICAL DISTANCES
    indCndDist = reshape(randsample(1:length(meanFocstmOptDstUnq),nCnd*1,1),[nCnd 1]);
    % OPTICAL DISTANCE DURING 1ST HALF OF TRIALS
    meanFocstmOptDst = meanFocstmOptDstUnq(indCndDist);
    % OPTICAL DISTANCE DURING 2ND HALF OF TRIALS
    meanFocstmOptDst2 = [];
    for i = 1:length(meanFocstmOptDst)
        meanFocstmOptDstTmp = meanFocstmOptDstUnq(meanFocstmOptDstUnq~=meanFocstmOptDst(i));
        meanFocstmOptDst2(i,:) = randsample(meanFocstmOptDstTmp,1);
    end
    % focStmOptDst = stepVals(indCnd(:,7)).*ones([nCnd 1]);
    % CURRENTLY CODE WANTS TO KNOW STEP CHANGE
    focStmOptDst = meanFocstmOptDst2-meanFocstmOptDst;
    AFCv = 1:length(focStmOptDst); % TRIAL ORDER (ALREADY RANDOMIZED)
end

if ~exist('sr')
   sr = [0 0]; 
end

% ----TURN IT OFF FOR NOW! ---
tcaCorrect = 0;
display('ARC_AFC: WARNING! For now, tcaCorrect has been manually set at 0');
% ----------------------------

if tcaCorrect==1
   load JnJ\TCAflsL; TCAfnmL=TCAfls{sn-1000,vs};
   load JnJ\TCAflsR; TCAfnmR=TCAfls{sn-1000,vs};
end

 if tcaCorrect==0; 
     LfarPower=opto(name_map('r_t_far')).control.getFocalPower.focal_power; 
     zL0=fnz0(LfarPower, double(ACL~=0)); zL=zL0(:,1:2);
 elseif tcaCorrect==1; 
     load(TCAfnmL, 'TCAp');
     zL=TCAp.sbjTCA; 
 end %zL=[0 0; 10 10; 20 20]; end; %TCAfnmL=TCAfnm; 
         
if tcaCorrect==0; 
    RfarPower=opto(name_map('r_t_far')).control.getFocalPower.focal_power; 
    zR0=fnz0(RfarPower, double(ACL~=0)); zR=zR0(:,3:4);
elseif tcaCorrect==1;
    load(TCAfnmR, 'TCAp'); 
    zR=TCAp.sbjTCA; 
end %zR=[0 0; 10 10; 20 20]; end; %TCAfnmR=TCAfnm; 
          
% if exist('AFCim0')==0 | isempty(AFCim0)==1; AFCfnm0='G:\My Drive\exp_bvams\code_repo\AFCim220510.mat'; load(AFCfnm0);  end %E optotype base 3 17secs to load


fprintf('Best shpere refraction: L = %f  , R = %f\n', sr(1), sr(2));
fprintf('Display power: L = %f  , R = %f\n',opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power);
fprintf('Trombone power: L = %f  , R = %f\n', opto(name_map('l_t_far')).control.getFocalPower.focal_power, opto(name_map('r_t_far')).control.getFocalPower.focal_power);
fprintf('Trombone position: L = %f  , R = %f\n',  zaber(name_map('l_trombone')).control.getposition,  zaber(name_map('r_trombone')).control.getposition);

disp(['Subject #' n2s(sn) ey ' EYE ACL' n2s(ACL) ' EXPERIMENT ' ex ' PRESS CTRL+C TO ABORT OR ANY OTHER KEY TO START']); KbWait([], 2);

[window1, window2, vbl0]=strt_psych0(screenNumber-2, screenNumber-1, 0);
        
cf=ones(3,2);
if ey(1)=='R'; cf(:,2)=0; elseif ey(1)=='L'; cf(:,1)=0; end

v0=focStmOptDst; 
v00=v0(AFCv,:); 
meanv00 = meanFocstmOptDst(AFCv);
rgb100 = rgb1(AFCv,:);
rgb200 = rgb2(AFCv,:);
% [im2L0, im2L1, im2R0, im2R1] = AFCtcaStmImg(AFCim0, AFCim1, zL, zR);
% im2 = 255.*insertText(zeros([500 500]),[125 75],'b','TextColor','green','BoxColor','black','FontSize',200);
% im2 = AFCwordStim('car',[500 500],[70 70; 160 70; 260 70]);
% im2(im2>0) = 255;
% im2 = flipud(im2);
AFCfls0=[filePath 'S' num2str(sn) 'V' num2str(vs) '_AFC_' ey 'ACL' n2s(ACL) '_' tme];

% imB = imread('G:\My Drive\exp_bvams\code_repo\imgs\vrn10_B_sd1.png');
imPattern = {};
for i = 1:6
    im1 = imread(['H:\Shared drives\CIVO_BVAMS\stimuli\word_image_0' num2str(i) '.png']);
    im1(im1>0) = 255;
    im1 = flipud(im1);   
    imPatternTmp = squeeze(im1(:,:,3));
    imPatternTmp = [zeros([30 size(imPatternTmp,2)]); imPatternTmp; zeros([30 size(imPatternTmp,2)])];
    imPatternTmp = [zeros([size(imPatternTmp,1) 30]) imPatternTmp zeros([size(imPatternTmp,1) 30])];
    imPattern{i} = imPatternTmp;
end

AFCp=AFCrgb(imPattern,rgb100,rgb200,v00, meanv00, sr, window1, window2, vs);    

AFCp.v0=v0;
AFCp.v00=v00;
AFCp.AFCv=AFCv;
AFCp.meanv00 = meanv00;
% AFCp.imfnm=AFCfnm0;
AFCp.rgb100 = rgb100;
AFCp.rgb200 = rgb200;

if sv == 1
%             AFCfls0=['JnJ\S' num2str(sn) 'V' num2str(vs) '_AFC_' ey 'ACL' n2s(ACL) '_' tme];
    save(AFCfls0, 'AFCp'); 
    load([filePath 'AFCfls' ey(1) '.mat'], 'AFCfls'); AFCfls{sn-1000,vs}=AFCfls0; 
    save([filePath 'AFCfls' ey(1) '.mat'], 'AFCfls'); 
end
         
opto(name_map('l_disp')).control.setFocalPower(14+sr(1));
opto(name_map('r_disp')).control.setFocalPower(14.3+sr(2));
% zaber(name_map('rotation')).move_deg(-3); %%-6400
[iLf iRf]=cwin3(imread('black.png'), imread('black.png') , cf, rc00, window2, window1);
clear LCAim;
sca;
