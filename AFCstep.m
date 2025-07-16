%%220508 AFC9f include TCA correction 
% function [v1 t3 dgs]=AFC8f(im_path, v0, sr, window1, window2)
function AFCp=AFCstep(im2L0, im2L1, im2R0, im2R1, v0, meanv0, sr, window1, window2)

%v0 vector of accomodations
%        control: [1×1 Zaber.AsciiDevice]
%     unit_scale: 2.1333e+03
%           move: @(x)zaber(name_map(ident)).control.moveabsolute(x)
%       move_deg: @(x)zaber(name_map(ident)).control.moveabsolute(x*zaber(name_map(ident)).unit_scale)
%        move_mm: []


% function [p1 p2 p3 p4 p5 p6 p7 p8 p9]=AFC0f(im_path, p0, window1, window2)
%p1(k0,:)=[power_dispL power_dispR powerL powerR rot];
global sz cf rc00 name_map zaber opto log

if size(v0,1)~=size(meanv0,1)
    error('AFCstep: v0 and meanv0 must have the same number of rows!');
end

nRmpSteps = 10; % number of steps for cosine ramp
nRmpSteps = round(nRmpSteps/2)*2; % make sure it's even number
tIntervalStm = 0.2; % how long to pause after each step
nFrmStmPlat = 10; % number of frames for which the stimulus plateaus

power_dispR=14.4+sr(2); %starting display power
power_dispL=14+sr(1); %starting display power
opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(1));
opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(1));

% wn=cwin0(img0, 'Stereo', cf, rc00, window1, window2);
[iLf0 iRf0]=cwin3(im2L0, im2R0, cf, rc00, window1, window2);

%tcpip
log.CRITICAL = 5;
log.ERROR = 4;
log.WARNING = 3;
log.INFO = 2;
log.DEBUG = 1;
log.LEVEL = log.DEBUG;
scene.enable_tcp=0;
scene.trial_num=1;

if scene.enable_tcp
    cmsg('TCP enabled', log.INFO, log.LEVEL);
    scene.tcp_socket = tcpip('169.229.228.75', 31000, 'NetworkRole', 'server');
    cmsg('Waiting for TCP socket connection...', log.INFO, log.LEVEL);
    fopen(scene.tcp_socket);
    cmsg('TCP connected!', log.INFO, log.LEVEL);
end

%[power_dispL power_dispR]=fcs_afc(window1, window2);
%tan=opp/adj
ipd=62./1000; %atand opposite/adjacent devi
%dgs=atand(1./(ipd.*(power_dispL+v0))); %degrees to rotate
%dgs=atand(ipd.*(power_dispL+v0)); %degrees to rotate
dgs=atand(ipd.*v0)./2-3; dgs(dgs<=-3)=-3; %degrees to rotate
dgs0=atand(ipd.*2)./2-3; dgs0(dgs0<=-3)=-3; %degrees to rotate

t0=zeros(length(v0), 6); t1=t0; t2=t0;
sinValuesAll = [];
% stage) 0stop 1record figure this out with Steve
disp('ready to start');  KbWait([], 2); 
for k0=1:size(v0,1)
      tSin = 0:(1/(nRmpSteps-1)):1; % support for sinusoidal modulation
      sinValues = [];
      for i = 1:size(v0,2)
         sinValuesTmp = (sin(2*pi*tSin-pi*1/2)+1).*0.5; % the modulation itself
         sinValuesTmp = meanv0(k0)+v0(k0,i).*[zeros([1 nFrmStmPlat]) sinValuesTmp(1:length(sinValuesTmp)/2) ones([1 nFrmStmPlat])];
         sinValues = [sinValues sinValuesTmp];
      end
      sinValuesAll(k0,:) = sinValues;
      %wn=cwin0(img1, 'Stereo', cf, rc00, window1, window2);
      [iLf1 iRf1]=cwin3(im2L1, im2R1, cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(k0));
      zaber(name_map('rotation')).move_deg(dgs(k0)); %%-6400

      %disp( n2s(v0(k0)));        
      fprintf('TRL= %f, L = %f  , R = %f , DEG = %f, Demand = %f\n' ,k0, opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power, (zaber(name_map('rotation')).control.getposition)./2.1333E3, v0(k0) );

      snd(250, 0.25); %pause(2.75);
      KbWait([], 2); 

      if scene.enable_tcp; send_tcp0(scene, 1); end; t0(k0,:)=clock;
      
      snd(1000, 0.2); pause(0.8);
      for i = 1:length(sinValues)
         opto(name_map('l_disp')).control.setFocalPower(power_dispL-sinValues(i));
         opto(name_map('r_disp')).control.setFocalPower(power_dispR-sinValues(i));
         pause(tIntervalStm);
      end
      snd(1000, 0.2); pause(0.8);
     
      if scene.enable_tcp; send_tcp0(scene, 0); end %stage) 0stop 1record
      %pause(3);
      %wn=cwin0(img0, 'Stereo', cf, rc00, window1, window2);
      [iLf0 iRf0]=cwin3(im2L0, im2R0, cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(k0));
    %           zaber(name_map('rotation')).move_deg(-3); %%-6400
      zaber(name_map('rotation')).move_deg(dgs0); %%-6400

      t1(k0,:)=clock;

      %snd(2000, 0.25);  pause(2.75);
      pause(3);

      %if scene.enable_tcp; send_tcp0(scene, 0); end %stage) 0stop 1record
      t2(k0,:)=clock;

    scene.trial_num=k0;
end

if scene.enable_tcp; fclose(scene.tcp_socket); end
AFCp.v1=power_dispR
AFCp.t3=cat(3, t0, t1, t2);
AFCp.dgs=dgs;
AFCp.imL0=im2L0;  AFCp.imR0=im2R0;
AFCp.imL1=im2L1;  AFCp.imR1=im2R1;
AFCp.sinValues = sinValuesAll;
