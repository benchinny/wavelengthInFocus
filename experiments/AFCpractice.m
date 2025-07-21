%%220508 AFC9f include TCA correction 
% function [v1 t3 dgs]=AFC8f(im_path, v0, sr, window1, window2)
function AFCp=AFCpractice(word1, v0, meanv0, sr, window1, window2,wordColor,recordAccommodation,bSizeCue)

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
    error('AFCpractice: v0 and meanv0 must have the same number of rows!');
end

if bSizeCue
    imSzPix = [320 320];
else
    imSzPix = [320 320];
end
if imSzPix(1)>320 || imSzPix(2)>320
    imSzPix = [320 320];
end
stimSizeScale = 3.3.*atand(1./(100./(meanv0(1).*0.8)))./atand(1./(100./(5))); % value of 3.3 gives a stimulus 3 degrees in size

nRmpSteps = 10; % number of steps for cosine ramp
nRmpSteps = round(nRmpSteps/2)*2; % make sure it's even number
tIntervalStm = 0.2; % how long to pause after each step
nFrmStmPlat = 10; % number of frames for which the stimulus plateaus

power_dispR=14.4+sr(2); %starting display power
power_dispL=14+sr(1); %starting display power
opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(1));
opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(1));

% MAKING WORD STIMULI
im0initial = AFCwordStimImproved(word1(1,:),imSzPix,wordColor);
im0initial(:,:,3) = circshift(squeeze(im0initial(:,:,3)),-15);

if bSizeCue
    im0 = [];
    im0(:,:,1) = imresize(squeeze(im0initial(:,:,1)),size(squeeze(im0initial(:,:,1))).*stimSizeScale);
    im0(:,:,2) = imresize(squeeze(im0initial(:,:,2)),size(squeeze(im0initial(:,:,2))).*stimSizeScale);
    im0(:,:,3) = imresize(squeeze(im0initial(:,:,3)),size(squeeze(im0initial(:,:,3))).*stimSizeScale);
else
    im0 = im0initial;
end

im0(im0>1) = 255;
im0 = flipud(im0);
im2L0 = im0; im2R0 = im0; 

% wn=cwin0(img0, 'Stereo', cf, rc00, window1, window2);
[iLf0 iRf0]=cwin3(im2L0, im2R0, cf, rc00, window1, window2);

%tcpip
log.CRITICAL = 5;
log.ERROR = 4;
log.WARNING = 3;
log.INFO = 2;
log.DEBUG = 1;
log.LEVEL = log.DEBUG;
scene.enable_tcp=recordAccommodation;
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
      stimSizeScale = 3.3.*atand(1./(100./(meanv0(k0).*0.8)))./atand(1./(100./(5)));
      % MAKING WORD STIMULI
      im0initial = AFCwordStimImproved(word1(k0,:),imSzPix,wordColor);
      im0initial(:,:,3) = circshift(squeeze(im0initial(:,:,3)),-15);
      if bSizeCue
          im0 = [];
          im0(:,:,1) = imresize(squeeze(im0initial(:,:,1)),size(squeeze(im0initial(:,:,1))).*stimSizeScale);
          im0(:,:,2) = imresize(squeeze(im0initial(:,:,2)),size(squeeze(im0initial(:,:,2))).*stimSizeScale);
          im0(:,:,3) = imresize(squeeze(im0initial(:,:,3)),size(squeeze(im0initial(:,:,3))).*stimSizeScale);
      else
          im0 = im0initial; 
      end
      im0(im0>1) = 255;
      im0 = flipud(im0);
      im2L0 = im0; im2R0 = im0;    
      % MAKING CONTINUOUS ACCOMMODATIVE STIMULUS
      tSin = 0:(1/(nRmpSteps-1)):1; % support for sinusoidal modulation
      sinValues = [];
      for i = 1:size(v0,2)
         sinValuesTmp = (sin(2*pi*tSin-pi*1/2)+1).*0.5; % the modulation itself
         sinValuesTmp = meanv0(k0)+v0(k0,i).*[zeros([1 nFrmStmPlat]) sinValuesTmp(1:length(sinValuesTmp)/2) ones([1 nFrmStmPlat])];
         sinValues = [sinValues imresize([meanv0 meanv0+v0(k0,i) meanv0+v0(k0,i)],[1 length(sinValuesTmp)],'nearest')];
      end
      sinValuesAll(k0,:) = sinValues;
      %wn=cwin0(img1, 'Stereo', cf, rc00, window1, window2);
      [iLf0 iRf0]=cwin3(im2L0, im2R0, cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(k0));
%      zaber(name_map('rotation')).move_deg(dgs(k0)); %%-6400

      %disp( n2s(v0(k0)));        
      fprintf('TRL= %f, L = %f  , R = %f , DEG = %f, Demand = %f\n' ,k0, opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power, (zaber(name_map('rotation')).control.getposition)./2.1333E3, v0(k0) );

      snd(250, 0.25); %pause(2.75);
  
      KbName('UnifyKeyNames');
 %     KbWait([], 2); 
      exitLoop = 0;
      % Control loop
      ListenChar(2);
      try
          opt_chk=0;
          while opt_chk==0
              [ keyIsDown, keyTime, keyCode ] = KbCheck;
              if keyIsDown
                  if keyCode(KbName('5'))
                      opt_chk = 1;
                      rsp = 5;
                      disp(['Subject response: ' num2str(rsp) ', right answer = ' word1(k0,:)]);
                      %end
                  elseif keyCode(KbName('4')) %| keyCode(KbName('Return'))
                      opt_chk = 1;
                      rsp = 4;
                      disp(['Subject response: ' num2str(rsp) ', right answer = ' word1(k0,:)]);
                  elseif keyCode(KbName('6')) %| keyCode(KbName('Return'))
                      opt_chk = 1;
                      rsp = 6;
                      disp(['Subject response: ' num2str(rsp) ', right answer = ' word1(k0,:)]);
                  elseif keyCode(KbName('Return')) %| keyCode(KbName('Return'))
                      exitLoop = 1;
                      opt_chk = 1;
                  else
                      disp('WRONG KEY'); snd(100, 0.25);
                  end     
              end
              % Key debounce routine, which waits for key to be released
              while keyIsDown
                  [ keyIsDown, keyTime, keyCode ] = KbCheck;
              end
          end
      catch ERROR
          if enable_optotunes
              for p = 1:6
                  opto(p).control.Close();
              end
          end
          if enable_trombones
              fclose(port);
              delete(port);
          end
          a = instrfind();
          if ~isempty(a) %isempty(Obj) returns logical 1 (true) if Obj is an empty ExptData object. Otherwise, it returns logical 0 (false). An empty ExptData object contains no data elements.
              fclose(a);
              delete(a)
              clear a
          end
          rethrow(ERROR)
      end
      ListenChar(0);
%      KbWait([], 2);      
      if exitLoop
         break; 
      end 

      if scene.enable_tcp; send_tcp0(scene, 1); end; t0(k0,:)=clock;
      
      snd(1000, 0.2); pause(2);
      
      snd(1000, 0.1); pause(0.2);
      snd(1000, 0.1); pause(0.8);
     
      if scene.enable_tcp; send_tcp0(scene, 0); end %stage) 0stop 1record
      %pause(3);
      %wn=cwin0(img0, 'Stereo', cf, rc00, window1, window2);
    %           zaber(name_map('rotation')).move_deg(-3); %%-6400
%      zaber(name_map('rotation')).move_deg(dgs0); %%-6400

      t1(k0,:)=clock;

      %snd(2000, 0.25);  pause(2.75);
%      pause(3);

      %if scene.enable_tcp; send_tcp0(scene, 0); end %stage) 0stop 1record
      t2(k0,:)=clock;

    scene.trial_num=k0;
end

if scene.enable_tcp; fclose(scene.tcp_socket); end
AFCp.v1=power_dispR
AFCp.t3=cat(3, t0, t1, t2);
AFCp.dgs=dgs;
AFCp.imL0=im2L0;  AFCp.imR0=im2R0;
AFCp.sinValues = sinValuesAll;
