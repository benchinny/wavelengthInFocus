%%220508 JnJ_AFC AFC9f include TCA correction. 
ex='RGB'; %ey=input('which eye? Right/Binc');
ey=input('which eye are we testing (type 1 for Right or 2 for Binocular)? ');
if ey==1; ey='Right'; elseif ey==2; ey='Binc'; end
filePath = 'H:\Shared drives\CIVO_BVAMS\data\ARC\';
vsIncrement1 = input(['Increment visit number? The current visit number is ' num2str(vs) ' (1 for yes, 0 for no)']);
vsIncrement2 = input(['Asking again: increment visit number? The current visit number is ' num2str(vs) ' (1 for yes, 0 for no)']);
if vsIncrement1==vsIncrement2
    vsIncrement = vsIncrement2;
else
    error('Start over: did not give same answer regarding whether to increment');
end

if exist('colorGroup','var')
    colorGroup1 = input(['Color group? The current group is ' colorGroupString '. Possibilities: 1, 2, or 3. ']);
    colorGroup2 = input(['Asking again: color group? The current group is ' num2str(colorGroupString) '. Possibilities: 1, 2, or 3. ']);
    if colorGroup1==colorGroup2
       colorGroup = colorGroup2;
       colorGroupString = [colorGroupString ', ' num2str(colorGroup(end))];
    else
       error('Start over: inconsistent answers about color group');
    end
else
    colorGroup1 = input(['Color group? The options are 1, 2, or 3. ']);
    colorGroup2 = input(['Asking again: color group? The options are 1, 2, or 3. ']);
    if colorGroup1==colorGroup2
       colorGroup = colorGroup2;
       colorGroupString = num2str(colorGroup);
    else
       error('Start over: inconsistent answers about color group');
    end
end

if vsIncrement>=1
   vs = vs+1;
end

meanFocstmOptDstUnq = [1.5 2.5 3.5]';
meanFocstmOptDstUnq = meanFocstmOptDstUnq*1.2255;

rng(8);
colorsAll = [0.569 0.000 1.000; ...
             0.472 0.000 0.815; ...
             0.569 0.432 1.000; ...
             0.569 0.334 1.000; ...
             0.432 0.334 1.000; ...
             0.327 0.334 1.000; ...
             0.569 0.334 0.740; ...   
             0.569 0.334 0.547; ...           
             0.432 0.000 1.000; ...
             0.327 0.000 1.000; ...
             0.569 0.000 0.740; ...   
             0.569 0.000 0.547; ...                  
            ];

colorScrambleInd = randperm(size(colorsAll,1));
colorsAllScrambled = colorsAll(colorScrambleInd,:);

if colorGroup==1
    nRepeats = 3;
    rgb1 = colorsAllScrambled(1:4,:);
elseif colorGroup==2
    nRepeats = 3;
    rgb1 = colorsAllScrambled(5:8,:);
elseif colorGroup==3
    nRepeats = 3;
    rgb1 = colorsAllScrambled(9:12,:);
end
       
% TRIALS
meanFocstmOptDst = repmat(imresize(meanFocstmOptDstUnq,[nRepeats*length(meanFocstmOptDstUnq) 1],'nearest'),[size(rgb1,1) 1]);
rgb1 = imresize(rgb1,[size(rgb1,1)*nRepeats*length(meanFocstmOptDstUnq) size(rgb1,2)],'nearest');
% AFCv IS JUST A VECTOR OF INDICES FOR RANDOMIZING TRIAL ORDER
AFCv = [];
for i = 1:size(unique(rgb1,'rows'),1)
    AFCv = [AFCv; length(meanFocstmOptDstUnq)*nRepeats*(i-1)+randsample(1:length(meanFocstmOptDstUnq)*nRepeats,length(meanFocstmOptDstUnq)*nRepeats)'];
end

if ~exist('sr')
   sr = [0 0]; 
end

% ----TURN OFF TCA CORRECTION ---
tcaCorrect = 0;

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

fprintf('Best sphere refraction: L = %f  , R = %f\n', sr(1), sr(2));
fprintf('Display power: L = %f  , R = %f\n',opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power);
fprintf('Trombone power: L = %f  , R = %f\n', opto(name_map('l_t_far')).control.getFocalPower.focal_power, opto(name_map('r_t_far')).control.getFocalPower.focal_power);
fprintf('Trombone position: L = %f  , R = %f\n',  zaber(name_map('l_trombone')).control.getposition,  zaber(name_map('r_trombone')).control.getposition);

disp(['Subject #' n2s(sn) ey ' EYE ACL' n2s(ACL) ' EXPERIMENT ' ex ' PRESS CTRL+C TO ABORT OR ANY OTHER KEY TO START']); KbWait([], 2);

[window1, window2, vbl0]=strt_psych0(screenNumber-2, screenNumber-1, 0);
        
cf=ones(3,2);
if ey(1)=='R'; cf(:,2)=0; elseif ey(1)=='L'; cf(:,1)=0; end

meanv00 = meanFocstmOptDst(AFCv);
rgb100 = rgb1(AFCv,:);

AFCfls0=[filePath 'S' num2str(sn) 'V' num2str(vs) '_AFC_' ey 'ACL' n2s(ACL) '_' tme];

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

AFCp=AFCrgbMain(imPattern,rgb100, meanv00, sr, window1, window2, vs);    

AFCp.AFCv=AFCv;
AFCp.meanv00 = meanv00;
% AFCp.imfnm=AFCfnm0;
AFCp.rgb100 = rgb100;

if sv == 1
    save(AFCfls0, 'AFCp'); 
    load([filePath 'AFCfls' ey(1) '.mat'], 'AFCfls'); AFCfls{sn-1000,vs}=AFCfls0; 
    save([filePath 'AFCfls' ey(1) '.mat'], 'AFCfls'); 
end

opto(name_map('l_disp')).control.setFocalPower(14+sr(1));
opto(name_map('r_disp')).control.setFocalPower(14.3+sr(2));
[iLf iRf]=cwin3(imread('black.png'), imread('black.png') , cf, rc00, window2, window1);
clear LCAim;
sca;
