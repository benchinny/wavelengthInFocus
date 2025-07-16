% BASIC INPUTS--SUBJECT AND VISIT NUMBERS
sn=input('Enter subject number?'); sn=sn+1000;

load('H:\Shared drives\CIVO_BVAMS\data\ARC\AFCflsR.mat');
vsEmpty = [];
if sn-1000 <= size(AFCfls,1)
    for i = 1:size(AFCfls,2)
        if isempty(AFCfls{sn-1000,i})
            vsEmpty = [vsEmpty i];
        end
    end
    vs=input(['Enter visit number? Next unused visit number is ' num2str(vsEmpty(1)) '\n']); 
else
    vsEmpty = 1;
    vs=input(['Enter visit number? Next unused visit number is ' num2str(vsEmpty(1)) '\n']); 
end

clear AFCfls;

filePath = 'H:\Shared drives\CIVO_BVAMS\data\ARC\';

tcaCorrect=0; % correct for TCA?
ACL=0; % is there an ACL lens in the system?
bFlexTrombones = true; % do we want to 'stretch' the trombones a bit before starting the experiment?

if ACL==1
   load([filePath 'LCAflsL']); LeftLCA=LCAfls{sn-1000,2};
   load([filePath 'LCAflsR']); RightLCA=LCAfls{sn-1000,2};

   load(LeftLCA,  'LCAp'); aFull=LCAp.Full; aHalf=LCAp.Half; sr=LCAp.sr; clear LCAp
else
    sr = [0 0]; % use this variable to adjust default optical distances for myopic or hyperopic subjects
end

power0=8; % this number is never actually used in any experiments
%Initialize power of left display optotune
power_dispL=14+sr(1); 
if ACL==0 powerL=opL(1); elseif ACL==1; powerL=aFull.p9(5,1); elseif ACL==2; powerL=aHalf.p9(5,1); end
deg=-3; % used to move the angle of the left exit pupil--not using it so far
trombone_power_left=powerL; % power both left trombone optotunes are set to

if ACL==1
   load(RightLCA, 'LCAp'); aFull=LCAp.Full; aHalf=LCAp.Half; sr=LCAp.sr; clear LCAp
else
    sr = [0 0]; % optical distance correction for myopic or hyperopic subjects
end

%Initialize power of right display optotune
power_dispR=14.3+sr(2); %aFull.p9(2,2) ;
if ACL==0; powerR=opR(1); elseif ACL==1; powerR=aFull.p9(5,2); elseif ACL==2; powerR=aHalf.p9(5,2); end
trombone_power_right=powerR; % power both right trombone optotunes are set to
     
% SET POWERS OF LEFT DISPLAY OPTOTUNES
opto(name_map('l_disp')).control.setFocalPower(power_dispL);
opto(name_map('l_disp')).control.getFocalPower.focal_power
                   
%CHANGE LEFT TROMBONE
% powerL=a18(5,1);
opto(name_map('l_t_near')).control.setFocalPower(powerL); %[6.8400 24.9500]
opto(name_map('l_t_far')).control.setFocalPower(powerL); %[7.2400 27.2950]
%zaber(name_map('l_trombone')).move(l_trombone_f(powerL));%[8.619 13.8]
if bFlexTrombones
   zaber(name_map('l_trombone')).move(l_trombone_f(power0));%[8.619 13.8]
end
zaber(name_map('l_trombone')).control.waitforidle();

[opto(name_map('l_t_near')).control.getFocalPower.focal_power...
opto(name_map('l_t_far')).control.getFocalPower.focal_power...
zaber(name_map('l_trombone')).control.getposition]

%ROTATION LEFT
% deg=-3;            
zaber(name_map('rotation')).move_deg(deg); %%-6400            
zaber(name_map('rotation')).control.getposition  

%LEFT RETICLE
% trombone_power_left=powerL;

fprintf('Display LEFT %0.7f \n',power_dispL);
fprintf('Trombone LEFT %0.7f \n',trombone_power_left);
fprintf('\n');
%output final reticle position
%[PupCtr_LtX,PupCtr_LtY,PupCtr_RtX,PupCtr_RtY]=findPupilCenter(trombone_power_left, trombone_power_right); %(LeftTromPwr, RightTromPwr);
     %    fprintf('New Left Pupil Center X = %f  , Y = %f\n',PupCtr_LtX,PupCtr_LtY);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%CHANGE DISPLAY OPTOTUNE RIGHT
% power_dispR=a18(2,2) ; %[6.9100 28.0750]
opto(name_map('r_disp')).control.setFocalPower(power_dispR);
opto(name_map('r_disp')).control.getFocalPower.focal_power
        
%CHANGE RIGHT TROMBONE
% powerR=a18(5,2); % opR(1);
% powerR=opR(1);
opto(name_map('r_t_near')).control.setFocalPower(powerR); %[6.8550 25.2250]
opto(name_map('r_t_far')).control.setFocalPower(powerR); %[7.2400 27.2950]
% zaber(name_map('r_trombone')).move(r_trombone_f(powerR)); %[7 20]
if bFlexTrombones
   zaber(name_map('r_trombone')).move(r_trombone_f(power0)); %[7 20]
end
zaber(name_map('r_trombone')).control.waitforidle();


[opto(name_map('r_t_near')).control.getFocalPower.focal_power...
opto(name_map('r_t_far')).control.getFocalPower.focal_power...
zaber(name_map('r_trombone')).control.getposition];
        
[PupCtr_LtX,PupCtr_LtY,PupCtr_RtX,PupCtr_RtY]=findPupilCenter(powerL, powerR); %(LeftTromPwr, RightTromPwr);
        
[window1, window2, vbl0]=strt_psych0(screenNumber-2, screenNumber-1, 0);
cf=ones(3,2); [iLf iRf]=cwin3(imread('texture0_1080_newfill_malt.png'), imread('texture0_1080_newfill_malt.png') , cf, rc00, window2, window1);                
fprintf('\nLx = %f  , Ly = %f, Rx = %f  , Ry = %f\n',PupCtr_LtX, PupCtr_LtY, PupCtr_RtX, PupCtr_RtY);         
if ACL~=0
 disp(['ACL=' n2s(ACL) ' Put ACL lens in! TCA=' n2s(tcaCorrect)]); KbWait([], 2); 
else 
    disp(['ACL=' n2s(ACL) ' Remove ACL lens out! TCA=' n2s(tcaCorrect) ]); KbWait([], 2);
end
zaber(name_map('l_trombone')).move(l_trombone_f(powerL));%[8.619 13.8]
zaber(name_map('l_trombone')).control.waitforidle();        
zaber(name_map('r_trombone')).move(r_trombone_f(powerR));%[8.619 13.8]
zaber(name_map('r_trombone')).control.waitforidle();          
                  
%RIGHT RETICLE
% trombone_power_right=powerR;
fprintf('\nDisplay RIGHT %0.7f\n',power_dispR);
fprintf('Trombone RIGHT %0.7f\n',trombone_power_right);
fprintf('\n');
%output final reticle position
% [PupCtr_LtX,PupCtr_LtY,PupCtr_RtX,PupCtr_RtY]=findPupilCenter(trombone_power_left, trombone_power_right); %(LeftTromPwr, RightTromPwr);
         %fprintf('New Right Pupil Center X = %f  , Y = %f\n',PupCtr_RtX,PupCtr_RtY);     

[iLf iRf]=cwin3(imread('black.png'), imread('black.png') , cf, rc00, window2, window1);
sca
