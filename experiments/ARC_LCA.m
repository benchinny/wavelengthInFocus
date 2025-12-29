
%%220105 based on JJVC_LCA run LCA for JnJ subjects includes half LCA calculation
ex='LCA'; ACLinit=ACL; ACL=0;
sn=input('Enter subject number? '); sn=sn+1000;
vs=input('Enter visit number? '); 
sr=input('input best spherical refraction in vector format [Left Right]: ');
ey=input('which eye are we testing (type 0 for Left or 1 for Right)? ');
if ey==0; ey='Left'; elseif ey==1; ey='Right'; end
filePath = 'G:\My Drive\exp_bvams\code_repo\ARC\';

%   power_disp_min  initial_power_disp max_power_disp mag_min   mag  mag_max rotation reptitions
%a0=[8               13.5                  17             0.9    1     1.3     -3       6];
a0=[8               16                  17             0.9    1     1.4     -3       6];
[PupCtr_LtX,PupCtr_LtY,PupCtr_RtX,PupCtr_RtY]=findPupilCenter(opL(1), opR(1)); %(LeftTromPwr, RightTromPwr);

opto(name_map('r_disp')).control.setFocalPower(14.4+sr(1));
opto(name_map('r_disp')).control.getFocalPower.focal_power;

opto(name_map('l_disp')).control.setFocalPower(14+sr(2));
opto(name_map('l_disp')).control.getFocalPower.focal_power;

            zaber(name_map('r_trombone')).move(r_trombone_f(8)); %[7 20]
            zaber(name_map('r_trombone')).control.waitforidle();
            zaber(name_map('l_trombone')).move(l_trombone_f(8));%[8.619 13.8]
            zaber(name_map('l_trombone')).control.waitforidle();
            
[window1, window2, vbl0]=strt_psych0(screenNumber-2, screenNumber-1, 0);
 cf=ones(3,2); [iLf iRf]=cwin3(imread('texture0_1080_newfill_malt.png'), imread('texture0_1080_newfill_malt.png') , cf, rc00, window2, window1);                
%ListenChar(2);                    
fprintf('Best shpere refraction: L = %f  , R = %f\n', sr(1), sr(2));
fprintf('Display power: L = %f  , R = %f\n',opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power);
 fprintf('Trombone power: L = %f  , R = %f\n', opto(name_map('l_t_far')).control.getFocalPower.focal_power, opto(name_map('r_t_far')).control.getFocalPower.focal_power);
 fprintf('Trombone position: L = %f  , R = %f\n',  zaber(name_map('l_trombone')).control.getposition,  zaber(name_map('r_trombone')).control.getposition);
fprintf('Pupil Center LX = %f  , LY = %f,  RX = %f  , RY = %f\n', PupCtr_LtX, PupCtr_LtY, PupCtr_RtX, PupCtr_RtY);
%fprintf('Rotation: %f degrees\n',rot);

% [window1, window2, vbl0]=strt_psych0(screenNumber-1, screenNumber, 0);
%  cf=ones(3,2); [iLf iRf]=cwin3(imread('texture0_1080_newfill_malt.png'), imread('texture0_1080_newfill_malt.png') , cf, rc00, window2, window1);                
% wn=cwin0(imread('texture0_1080_newfill_malt.png'), 'Stereo', cf, rc00, window1, window2);
%disp(['Subject #' n2s(sn) ey ' EYE ACL' n2s(ACL) ' EXPERIMENT ' ex ' PRESS CTRL+C TO ABORT OR ANY OTHER KEY TO START']); KbWait([], 2);
if ACL~=0
    disp('ACL~=0!!!')
else
    disp(['Subject #' n2s(sn)  'Visit #' n2s(vs) ey ' EYE ACL' n2s(ACL) ' EXPERIMENT ' ex ' PRESS CTRL+C TO ABORT OR ANY OTHER KEY TO START']); %KbWait([], 2);
end
disp(['SYSTEM IS IN ACL CHANGE POSITION, PRESS ANY BUTTON FOR NEUTRAL POSITION AND EXPERIMENT START'])
KbWait([], 2); 
 %[window1, window2, vbl0]=strt_psych0(screenNumber-1, screenNumber, 0);
 %cf=ones(3,2); %[iLf iRf]=cwin3(imread('texture0_1080_newfill_malt.png'), imread('texture0_1080_newfill_malt.png') , cf, rc00, window2, window1);                

        %%input a output b
        if ey(1)=='R'
            cf(:,1)=0;
        elseif ey(1)=='L'
            cf(:,2)=0;
        end

        %initial power=13.5 AR
        
        LCAim{1}='vrn10_R_sd1.png'; %'text2.png'; %
        LCAim{2}='vrn10_G_sd1.png';
        LCAim{3}='vrn10_B_sd1.png'; %'opt1.png'; %
        %           LCAim{1}='E500by500inv.png'
        %           LCAim{2}='E500by500inv.png'
        %           LCAim{3}='E500by500inv.png'
        %[a10, a11, a12, a13, a14, a15, a16, a17, a18, a19]=LCA16f(LCAim, a0, window1, window2); %p1(a10):data, p2:avg, p3:sd, (a13)p4:sys+eye, p5(a14):EyeLCA, p6:trom pwr
        [aFull, aHalf]=LCA16f(LCAim, a0, window1, window2); %p1(a10):data, p2:avg, p3:sd, (a13)p4:sys+eye, p5(a14):EyeLCA, p6:trom pwr

        %a19 is HalfLCA parameters!
%         [PupCtr_LtX,PupCtr_LtY,PupCtr_RtX,PupCtr_RtY]=findPupilCenter(a18(5,1), a18(5,2)); %(LeftTromPwr, RightTromPwr);
%         fprintf('New Right Pupil Center X = %f  , Y = %f\n',PupCtr_RtX,PupCtr_RtY);
%         fprintf('New Left Pupil Center X = %f  , Y = %f\n',PupCtr_LtX,PupCtr_LtY);
        
        
        %findz0;
        
        %z0=fnz0(a18(5,:), ACL); %find system TCA
        
        VRCMprm; %save varichrome parameters into VRCMp
        LCAp=VRCMp; LCAp.Full=aFull; LCAp.Half=aHalf; LCAp.im=LCAim; LCAp.sr=sr;
        %if sv==1; save(['data\S' num2str(sn) ey '_' ex '_ACL' n2s(ACL) '_' tme], 'a0', 'a10', 'a11', 'a12', 'a13', 'a14', 'a15', 'a16', 'a17', 'a18', 'a19' ,'LCAim','PupCtr_RtX','PupCtr_RtY','PupCtr_LtX','PupCtr_LtY','z0','ACL', 'LCAp'); end
        if sv==1; 
            LCAfls0=[filePath 'S' num2str(sn) 'V' num2str(vs) ey '_' ex '_ACL' n2s(ACL) '_' tme];
            save(LCAfls0, 'LCAp'); 
            load([filePath 'LCAfls' ey(1) '.mat'], 'LCAfls'); LCAfls{sn-1000,vs}=LCAfls0; 
            save([filePath 'LCAfls' ey(1) '.mat'], 'LCAfls'); 
        end

%         save('JnJ\SPTtmp', 'aFull', 'aHalf', '-append');
                save([filePath 'SPTtmp.mat'], 'LCAp');

        %input(['\n PLEASE UPDATE PUPIL CENTER IN LABVIEW \n' 'rtX: ' n2s2(PupCtr_RtX) ' rtY: ' n2s2(PupCtr_RtY) ' ltX: ' n2s2(PupCtr_LtX) ' ltY: ' n2s2(PupCtr_LtY)])
        
        %a18
        opto(name_map('l_disp')).control.setFocalPower(14+sr(1));
        opto(name_map('r_disp')).control.setFocalPower(14.4+sr(2));
        zaber(name_map('rotation')).move_deg(-3); %%-6400
        [iLf iRf]=cwin3(imread('black.png'), imread('black.png') , cf, rc00, window2, window1);
        ACL=ACLinit;
        sca
        
        if ey(1)=='R'
           
        fprintf('Chromatic Difference of Refraction :\n Red-Blue = %4.2f  ,  Green-Blue = %4.2f  ,Red-Green = %4.2f\n', LCAp.Full.p5(1,2), LCAp.Full.p5(2,2), LCAp.Full.p5(3,2));

        fprintf('Mean Display Focus +-STD :\n Red = %4.2f +- %4.2f  ,  Green = %4.2f +- %4.2f  ,Blue = %4.2f +- %4.2f\n', LCAp.Full.p2(1,2), LCAp.Full.p3(1,2), LCAp.Full.p2(2,2), LCAp.Full.p3(2,2), LCAp.Full.p2(3,2), LCAp.Full.p3(3,2));

        elseif ey(1)=='L'
%         
        fprintf('Chromatic Difference of Refraction :\n Red-Blue = %4.2f  ,  Green-Blue = %4.2f  ,Red-Green = %4.2f\n', LCAp.Full.p5(1,1), LCAp.Full.p5(2,1), LCAp.Full.p5(3,1));

        fprintf('Mean Display Focus +-STD :\n Red = %4.2f +- %4.2f  ,  Green = %4.2f +- %4.2f  ,Blue = %4.2f +- %4.2f\n', LCAp.Full.p2(1,1), LCAp.Full.p3(1,1), LCAp.Full.p2(2,1), LCAp.Full.p3(2,1), LCAp.Full.p2(3,1), LCAp.Full.p3(3,1));
        end
% %%210604 JJVC_LCA run LCA task
% ex='LCA';
% %   power_disp_min  initial_power_disp max_power_disp mag_min   mag  mag_max rotation reptitions
% a0=[8               13.5                  17             0.9    1     1.3     -3       6];
% a0=[8               13.5                  17             0.9    1     1.4     -3       6];
% 
% 
% fprintf('Display power: L = %f  , R = %f\n',opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power);
%  fprintf('Trombone power: L = %f  , R = %f\n', opto(name_map('l_t_far')).control.getFocalPower.focal_power, opto(name_map('r_t_far')).control.getFocalPower.focal_power);
%  fprintf('Trombone position: L = %f  , R = %f\n',  zaber(name_map('l_trombone')).control.getposition,  zaber(name_map('r_trombone')).control.getposition);
% fprintf('Pupil Center LX = %f  , LY = %f,  RX = %f  , RY = %f\n', PupCtr_LtX, PupCtr_LtY, PupCtr_RtX, PupCtr_RtY);
% %fprintf('Rotation: %f degrees\n', rot);
% 
% disp(['Subject #' n2s(sn) ey ' EYE ACL' n2s(ACL) ' EXPERIMENT ' ex ' PRESS CTRL+C TO ABORT OR ANY OTHER KEY TO START']); KbWait([], 2);
% 
% [window1, window2, vbl0]=strt_psych0(screenNumber-1, screenNumber, 0);
% 
%         %%input a output b
%         cf=ones(3,2);
% 
%         %initial power=13.5 AR
%         
%         LCAim{1}='vrn10_R_sd1.png'; %'text2.png'; %
%         LCAim{2}='vrn10_G_sd1.png';
%         LCAim{3}='vrn10_B_sd1.png'; %'opt1.png'; %
%         %           LCAim{1}='E500by500inv.png'
%         %           LCAim{2}='E500by500inv.png'
%         %           LCAim{3}='E500by500inv.png'
%         %[a10, a11, a12, a13, a14, a15, a16, a17, a18]=LCA12f(LCAim, a0, window1, window2); %p1(a10):data, p2:avg, p3:sd, (a13)p4:sys+eye, p5(a14):EyeLCA, p6:trom pwr
%         [a10, a11, a12, a13, a14, a15, a16, a17, a18]=LCA14f(LCAim, a0, window1, window2); %p1(a10):data, p2:avg, p3:sd, (a13)p4:sys+eye, p5(a14):EyeLCA, p6:trom pwr
% 
%         [PupCtr_LtX,PupCtr_LtY,PupCtr_RtX,PupCtr_RtY]=findPupilCenter(a18(5,1), a18(5,2)); %(LeftTromPwr, RightTromPwr);
%         fprintf('New Right Pupil Center X = %f  , Y = %f\n',PupCtr_RtX,PupCtr_RtY);
%         fprintf('New Left Pupil Center X = %f  , Y = %f\n',PupCtr_LtX,PupCtr_LtY);
%         
%         
%         %findz0;
%         
%         z0=fnz0(a18(5,:), ACL); %find system TCA
%         VRCMprm; %save varichrome parameters into VRCMp
%         LCAp=VRCMp;
%         if sv==1; save(['data\S' num2str(sn) ey '_' ex '_ACL' n2s(ACL) '_' tme], 'a0', 'a10', 'a12', 'a13', 'a14', 'a15', 'a16', 'a17', 'a18', 'a11','LCAim','PupCtr_RtX','PupCtr_RtY','PupCtr_LtX','PupCtr_LtY','z0','ACL', 'LCAp'); end
% 
% %         if sv==1; save(['data\S' num2str(sn) '_LCA_' ey 'ACL' n2s(ACL) '_' tme], 'a0', 'a10', 'a12', 'a13', 'a14', 'a15', 'a16', 'a17', 'a18', 'a11','LCAim','PupCtr_RtX','PupCtr_RtY','PupCtr_LtX','PupCtr_LtY','z0','ACL'); end
%         save('data\SPTtmp', 'a18', '-append');
%         
%         input(['\n PLEASE UPDATE PUPIL CENTER IN LABVIEW \n' 'rtX: ' n2s2(PupCtr_RtX) ' rtY: ' n2s2(PupCtr_RtY) ' ltX: ' n2s2(PupCtr_LtX) ' ltY: ' n2s2(PupCtr_LtY)])
%         
%         a18
%         sca
%         
        