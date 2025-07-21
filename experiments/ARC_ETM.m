  %%220105 JnJ Accommodation
ex='ETM'; %ey=input('which eye? Right/Binc');

ey=input('which eye are we testing (type 1 for Right or 2 for Binocular)? ');
if ey==1; ey='Right'; elseif ey==2; ey='Binc'; end

% TCAfnmL='S1004V4Left_TCA_ACL2_2204121430.mat';
% TCAfnmR='S1004V4Right_TCA_ACL2_2204121434.mat'; 
% 
if tcaCorrect==1
load JnJ\TCAflsL; TCAfnmL=TCAfls{sn-1000,vs};
load JnJ\TCAflsR; TCAfnmR=TCAfls{sn-1000,vs};
end
%        input e output f
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% ETM params
        %y0=zeros(size(y0)); d1=5; [r00 r11]=psf_ETM1(y0, y1,STK11(w2),d1); disp(r00); disp(r11);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %bxy=[-10.3 10.7; 10.5 -11]
        d1=5; % 5trials to drop
        %n0=4;%4; % number of experiments/repetitions
        n1=6; %40;% 40 number of trials
    cntrst=1; %0.2; %1;%cntrst=1
    %tcaCorrect=1;
    acV=[-0.5:0.5:3];
    %y0=acV([7	8	3	5	6	2	4	1]); %y0=y0([1 2])
    %y0=acV([2 6 8 1 7 5 3 4]); %y0=y0([1 2])
    y0=acV(ETMv); %y0=y0([1 2])

cf=ones(3,2); n0=4; cg='w';
if ey(1)=='R'; cf(:,2)=0; elseif ey(1)=='L'; cf(:,1)=0; end

      %    !!!!!!!!!!!!!! UPDATE TCA and LCA files!!!!!!!!!!!!!!
        w0=13;% initial size of the optotype, stroke size= 2^w, marty 12
        sz = [1080, 1920]; % size of screen
        d0=0.5; % 0.5length of optotype presentation in secs, AR 0.5sec, MB 2 sec
        % Get the screen numbers
        screens = Screen('Screens');
        screenNumber = max(screens);
        load cal_val; %cf=[RB./RR LB./LR]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
        %load TCA file
        %TCAfnm='S1001V2Left_TCA_ACL0_2201111534.mat'; load(TCAfnm, 'TCAp');
        %TCAfnm='S1001V1Left_TCA_ACL1_2201062350.mat'; load(TCAfnm, 'TCAp');
         %load(TCAfnmL, 'TCAp');
         %if tcaCorrect==0; zL=[0 0; -0.2665 0.5690; -1.2130 -1.3906]; elseif tcaCorrect==1; load(TCAfnmL, 'TCAp'); zL=TCAp.sbjTCA; end; %TCAfnmL=TCAfnm; 
         if tcaCorrect==0; 
             LfarPower=opto(name_map('r_t_far')).control.getFocalPower.focal_power; 
             zL0=fnz0(LfarPower, double(ACL~=0)); zL=zL0(:,1:2);
         elseif tcaCorrect==1; load(TCAfnmL, 'TCAp'); zL=TCAp.sbjTCA; end; %TCAfnmL=TCAfnm; 

%          load('.mat', 'TCAp');
%         if tcaCorrect==0; zL=TCAp.sysTCA; elseif tcaCorrect==1; zL=TCAp.sbjTCA; end; %TCAfnmL=TCAfnm; 
%         z0L=[0 0; -10 -10; 10 10]; z1L=z0L;
%         z0L=[0 0; 0 0; 0 0]; z1L=z0L;

        %TCAfnm='S1001V2Right_TCA_ACL0_2201111531.mat'; load(TCAfnm, 'TCAp');
        %TCAfnm='S1001V1Right_TCA_ACL1_2201070001.mat'; load(TCAfnm, 'TCAp');
         %load(TCAfnmR, 'TCAp');
%                 if tcaCorrect==0; zR=TCAp.sysTCA; elseif tcaCorrect==1; zR=TCAp.sbjTCA; end; %TCAfnmR=TCAfnm; 
          if tcaCorrect==0; 
                RfarPower=opto(name_map('r_t_far')).control.getFocalPower.focal_power; 
                zR0=fnz0(RfarPower, double(ACL~=0)); zR=zR0(:,3:4);
          elseif tcaCorrect==1; load(TCAfnmR, 'TCAp'); zR=TCAp.sbjTCA; end; %TCAfnmR=TCAfnm; 

%         z0R=[0 0; 0 0; 0 0]; z1R=z0R;
%          z0R=[0 0; -10 -10; 10 10]; z1R=z0R;
        %%%%%%%%%%%%%CHECK RETICLE%%%%%%%%%%%%%%
        trombone_power_left=opto(name_map('l_t_far')).control.getFocalPower.focal_power; %a18(5,1);
        trombone_power_right=opto(name_map('r_t_far')).control.getFocalPower.focal_power; %a18(5,2);
        fprintf('Best shpere refraction: L = %f  , R = %f\n', sr(1), sr(2));
        fprintf('Display LEFT %0.7f \n',opto(name_map('l_disp')).control.getFocalPower.focal_power);
        fprintf('Trombone LEFT %0.7f \n',trombone_power_left);
        fprintf('\nDisplay RIGHT %0.7f\n',opto(name_map('r_disp')).control.getFocalPower.focal_power);
        fprintf('Trombone RIGHT %0.7f\n',trombone_power_right);
        fprintf('\n');
        %output final reticle position
        [PupCtr_LtX,PupCtr_LtY,PupCtr_RtX,PupCtr_RtY]=findPupilCenter(trombone_power_left, trombone_power_right); %(LeftTromPwr, RightTromPwr);
        fprintf('New Right Pupil Center X = %f  , Y = %f\n',PupCtr_RtX,PupCtr_RtY);
        fprintf('New Left Pupil Center X = %f  , Y = %f\n',PupCtr_LtX,PupCtr_LtY);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %NOTE: NEED TO SELECT APPROPRAITE VALUE OF ACL for findz0 to work
        %correctly!
        %          input(['\n WAS ACL VALUE UPDATED CORRECTLY?? Currently it is: ' n2s2(ACL)]);

        %if exist('E10c')==0 | isempty(E10c)==1; ETMfnm0='E_b2s5_p2_p14_o4_ms10_210320.mat'; load(ETMfnm0);  end %E optotype base 3 17secs to load
        if exist('E10c')==0 | isempty(E10c)==1; ETMfnm0='E_b2s5_p5_p17_o4_ms10_210621.mat'; load(ETMfnm0);  end %E optotype base 3 17secs to load

                    %if exist('E10c')==0 | isempty(E10c)==1;
        %ETMfnm0='E_b2s5_p10_p22_o4_ms10_210608.mat';load(ETMfnm0);end %JOL

        % 'E_b2s5_p2_p14_o4_ms10_210320.mat' acuities: 6 to  32: namely 6.0803    6.9844    8.0230    9.2160   10.5864   12.1606   13.9688 16.0460   18.4320   21.1728   24.3212   27.9377   32.0920
        % 'E_b2s2_p1_p11_o4_ms10_210319.mat' acuities:6 to 417
        % E_b2s10_p1_p13_o4_ms10_201207.mat buggy
        % 'E_b2s6_p2_p14_o4_ms10_210318.mat' 'E_b2s10_p1_p13_o4_ms10_210316.mat'
        disp(['Subject #' n2s(sn) ey ' EYE ACL' n2s(ACL) ' EXPERIMENT ' ex ' PRESS CTRL+C TO ABORT OR ANY OTHER KEY TO START']); KbWait([], 2);
        [window1, window2, vbl0]=strt_psych0(screenNumber-2, screenNumber-1, 1); %start psychtoolbox
        
        %presented in command window: trial #, base-power, stroke size, optotype size, time to generate optotype
        %if cg=='g'; z1L=zeros(size(z1L)); z0=z1; end; %zero TCA offset for green stim
        if cg=='g'; z1L=zeros(3,2); z0L=z1L; z0R=z0L; z1R=z0R; end %zero TCA offset for green stim

        %[y0 y1 w2 y5]=ETM16f([n0 n1 w0 d0], z1, z0, E10c, STK, window1, window2);
        %[y1 w2 y5]=ETM18f([n0 n1 w0 d0 tcaCorrect], y0, z1R, z0R, z1L, z0L, E10c, STK, sr, window1, window2);

        ETMp=ETMrgb([n0 n1 w0 d0], y0, zL, zR, E10c, STK, sr, window1, window2);
        ETMp.y0=y0;
        ETMp.ETMv=ETMv;
        ETMp.acV=acV;

        %y0=zeros(size(y0)); d1=5; [r00 r11]=psf_ETM(y0, y1, w2, d1, STK00);
        %%y5=save final varichrome settings
        c0=0.00384; %deg/pixels
        c1=1./(c0.*60); %pixels/minute
        w3=ETMp.w2; rc=size(ETMp.w2); for r=1:rc(1); for c=1:rc(2); w3(r, c)=STK(ETMp.w2(r,c)); end; end;
        w4=20.*w3./c1; 20.*STK./c1
        ETMp.STK11=STK00; 
        
        ETMprm;
        ETMp.opt_fnm=ETMfnm0;
%         ETMp.sysTcaL=z0L; ETMp.sysTcaR=z0R;
%         ETMp.sbjTcaL=z1L; ETMp.sbjTcaL=z1L;
        ETMp.tcaL=zL; ETMp.tcaR=zR; ETMp.tcaCorrect=tcaCorrect;
        if tcaCorrect==1
         ETMp.TCAfnmL=TCAfnmL; clear TCAp;
         ETMp.TCAfnmR=TCAfnmR;
        end
        %if sv == 1; save(['data\S' num2str(sn) '_ETM' cg '_' ey 'ACL' n2s(ACL) '_' tme], 'y0', 'y1', 'w2', 'z0', 'z1', 'ETMfnm0', 'y5', 'VRCMp'); end
        %if sv == 1; save(['data\S' num2str(sn) ey '_' ex cg '_ACL' n2s(ACL) '_C' n2s(round(100.*cntrst)) '_' tme], 'y0', 'y1', 'w2', 'y5','STK11', 'ETMp'); end
%         if sv == 1; save(['JnJ\S' num2str(sn) 'V' num2str(vs) ey '_' ex cg '_ACL' n2s(ACL) '_C' n2s(round(100.*cntrst)) '_' tme], 'y0', 'y1', 'w2', 'y5','STK11', 'ETMp'); end
       % if sv == 1; save(['JnJ\S' num2str(sn) 'V' num2str(vs) ey '_' ex cg '_ACL' n2s(ACL) '_C' n2s(round(100.*cntrst)) '_' tme], 'ETMp'); end
if sv == 1;
            ETMfls0=['JnJ\S' num2str(sn) 'V' num2str(vs) ey '_' ex cg '_ACL' n2s(ACL) '_C' n2s(round(100.*cntrst)) '_' tme];
    save([ETMfls0 '2'], 'ETMp'); 
    load(['JnJ\ETMfls' ey(1) '.mat'], 'ETMfls'); ETMfls{sn-1000,vs}=ETMfls0; save(['JnJ\ETMfls2' ey(1) '.mat'], 'ETMfls'); 
end
        
        
        
        %% z1 & z2 are TCA values used (mean & std)
          opto(name_map('l_disp')).control.setFocalPower(14+sr(1));
          opto(name_map('r_disp')).control.setFocalPower(14.4+sr(2));
          zaber(name_map('rotation')).move_deg(-3); %%-6400
        [iLf iRf]=cwin3(imread('black.png'), imread('black.png') , cf, rc00, window2, window1);
        sca
