%%220508 JnJ_AFC AFC9f include TCA correction. 
ex='ACU'; %ey=input('which eye? Right/Binc');
ey=1;
if ey==1; ey='Right'; elseif ey==2; ey='Binc'; end
filePath = 'G:\My Drive\exp_bvams\code_repo\ARC\';
vsIncrement = input(['Increment visit number? The current visit number is ' num2str(vs) ' (1 for yes, 0 for no)']);
vsIncrement = input(['Asking again: increment visit number? The current visit number is ' num2str(vs) ' (1 for yes, 0 for no)']);
% note: GIVE PEOPLE MORE TIME TO ACCOMMODATE TO THE NEW STIMULUS
if vsIncrement>=1
   vs = vs+1;
end

meanFocstmOptDst = [1.5 2.5 3.5]*1.2255;
meanFocstmOptDst = meanFocstmOptDst(randperm(length(meanFocstmOptDst)));
rgb  = [
        0.569 0.432 1.00; ...
        ];
indScrambleRgb = randperm(size(rgb,1));
rgb = rgb(indScrambleRgb,:);
stimSizePix = [125];
stimSizePix = stimSizePix(indScrambleRgb);

focStmOptDstIncr = [0];
focStmOptDstIncr = focStmOptDstIncr.*1.2255;
trlPerLvl = 20;

% DEFAULT NO TCA CORRECTION
LfarPower=opto(name_map('r_t_far')).control.getFocalPower.focal_power; 
zL0=fnz0(LfarPower, double(ACL~=0)); zL=zL0(:,1:2);
% DEFAULT NO TCA CORRECTION
RfarPower=opto(name_map('r_t_far')).control.getFocalPower.focal_power; 
zR0=fnz0(RfarPower, double(ACL~=0)); zR=zR0(:,3:4);

fprintf('Display power: L = %f  , R = %f\n',opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power);
fprintf('Trombone power: L = %f  , R = %f\n', opto(name_map('l_t_far')).control.getFocalPower.focal_power, opto(name_map('r_t_far')).control.getFocalPower.focal_power);
fprintf('Trombone position: L = %f  , R = %f\n',  zaber(name_map('l_trombone')).control.getposition,  zaber(name_map('r_trombone')).control.getposition);

disp(['Subject #' n2s(sn) ey ' EYE ACL' n2s(ACL) ' EXPERIMENT ' ex ' PRESS CTRL+C TO ABORT OR ANY OTHER KEY TO START']); KbWait([], 2);

[window1, window2, vbl0]=strt_psych0(screenNumber-2, screenNumber-1, 0);
        
cf=ones(3,2);
if ey(1)=='R'; cf(:,2)=0; elseif ey(1)=='L'; cf(:,1)=0; end

AFCfls0=[filePath 'S' num2str(sn) 'V' num2str(vs) '_AFC_' ey 'ACL' n2s(ACL) '_' tme];

imPattern = [];
im1 = AFCwordStimImproved('sun',[320 320],'green');
im1(im1>0) = 255;
im1 = flipud(im1);   
imPatternTmp = squeeze(im1(:,:,2));
imPatternTmp = circshift(imPatternTmp,15,1);
imPattern(:,:,1) = imresize(imPatternTmp,[480 480]);
im2 = AFCwordStimImproved('sea',[320 320],'green');
im2(im2>0) = 255;
im2 = flipud(im2); 
imPatternTmp = squeeze(im2(:,:,2));
imPatternTmp = circshift(imPatternTmp,15,1);
imPattern(:,:,2) = imresize(imPatternTmp,[480 480]);
im3 = AFCwordStimImproved('ace',[320 320],'green');
im3(im3>0) = 255;
im3 = flipud(im3); 
imPatternTmp = squeeze(im3(:,:,2));
imPatternTmp = circshift(imPatternTmp,15,1);
imPattern(:,:,3) = imresize(imPatternTmp,[480 480]);
im4 = AFCwordStimImproved('one',[320 320],'green');
im4(im4>0) = 255;
im4 = flipud(im4); 
imPatternTmp = squeeze(im4(:,:,2));
imPatternTmp = circshift(imPatternTmp,15,1);
imPattern(:,:,4) = imresize(imPatternTmp,[480 480]);
AFCp=ARCacuityAfterAccFuncBlock(imPattern,rgb,meanFocstmOptDst,focStmOptDstIncr, window1, window2, trlPerLvl, stimSizePix);    

if sv == 1
    save(AFCfls0, 'AFCp'); 
    load([filePath 'AFCfls' ey(1) '.mat'], 'AFCfls'); AFCfls{sn-1000,vs}=AFCfls0; 
    save([filePath 'AFCfls' ey(1) '.mat'], 'AFCfls'); 
end

opto(name_map('l_disp')).control.setFocalPower(14);
opto(name_map('r_disp')).control.setFocalPower(14.3);
% zaber(name_map('rotation')).move_deg(-3); %%-6400
[iLf iRf]=cwin3(imread('black.png'), imread('black.png') , cf, rc00, window2, window1);
clear LCAim;
sca;
