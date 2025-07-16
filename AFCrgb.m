%%220508 AFC9f include TCA correction 
% function [v1 t3 dgs]=AFC8f(im_path, v0, sr, window1, window2)
function AFCp=AFCrgb(imPattern, rgb0, rgb1, v0, meanv0, sr, window1, window2, vs)

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
    error('AFCrgb: v0 and meanv0 must have the same number of rows!');
end

nRmpSteps = 10; % number of steps for cosine ramp
nRmpSteps = round(nRmpSteps/2)*2; % make sure it's even number
tIntervalStm = 0.2; % how long to pause after each step
nFrmStmPlat = 10; % number of frames for which the stimulus plateaus

power_dispR=14.3+sr(2); %starting display power
power_dispL=14+sr(1); %starting display power
opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(1));
opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(1));

if length(imPattern)>1
   indImPattern = randsample(1:length(imPattern),1);
else
   indImPattern = 1;
end
im2R0(:,:,1) = imPattern{indImPattern}.*rgb0(1,1);
im2R0(:,:,2) = imPattern{indImPattern}.*rgb0(1,2);
im2R0(:,:,3) = imPattern{indImPattern}.*rgb0(1,3);

% wn=cwin0(img0, 'Stereo', cf, rc00, window1, window2);
[iLf0 iRf0]=cwin3(im2R0, im2R0, cf, rc00, window1, window2);

%tcpip
log.CRITICAL = 5;
log.ERROR = 4;
log.WARNING = 3;
log.INFO = 2;
log.DEBUG = 1;
log.LEVEL = log.DEBUG;
scene.enable_tcp=1;
scene.trial_num=1;

if scene.enable_tcp && ~ismember('tcp_socket', who('global'))
    cmsg('TCP enabled');
    global tcp_socket;
    tcp_socket = tcpserver('169.229.228.57',31000);
    cmsg('Waiting for TCP socket connection...');
    bConnected = false;
    while ~bConnected
        tcpStatus = tcp_socket.Connected;
        pause(0.1);
        if tcpStatus
            bConnected=true;
        end
    end
    % fopen(tcp_socket);
    cmsg('TCP connected!');    
else
    global tcp_socket;
end

%[power_dispL power_dispR]=fcs_afc(window1, window2);
%tan=opp/adj
ipd=62./1000; %atand opposite/adjacent devi
%dgs=atand(1./(ipd.*(power_dispL+v0))); %degrees to rotate
%dgs=atand(ipd.*(power_dispL+v0)); %degrees to rotate
dgs=atand(ipd.*v0)./2-3; dgs(dgs<=-3)=-3; %degrees to rotate
dgs0=atand(ipd.*2)./2-3; dgs0(dgs0<=-3)=-3; %degrees to rotate

t0=zeros(length(v0), 6); t1=t0; t2=t0; tChange1 = t0; tChange2 = t0; tRealEnd = t0;
sinValuesAll = [];
% stage) 0stop 1record figure this out with Steve
disp('ready to start');  KbWait([], 2); 
indImPatternAll = [];
for k0=1:size(v0,1)
      tSin = 0:(1/(nRmpSteps-1)):1; % support for sinusoidal modulation
      sinValues = [];
      for i = 1:size(v0,2)
         sinValuesTmp = (sin(2*pi*tSin-pi*1/2)+1).*0.5; % the modulation itself
         sinValuesTmp = meanv0(k0)+v0(k0,i).*[zeros([1 nFrmStmPlat]) sinValuesTmp(1:length(sinValuesTmp)/2) ones([1 nFrmStmPlat])];
         sinValues = [sinValues imresize([meanv0(k0) meanv0(k0)+v0(k0,i)],[1 length(sinValuesTmp)],'nearest')];
      end
      sinValuesAll(k0,:) = sinValues;
      if length(imPattern)>1
        indImPattern = randsample(1:length(imPattern),1);
      else
        indImPattern = 1;
      end
      indImPatternAll(end+1,:) = indImPattern;
      im2R0(:,:,1) = imPattern{indImPattern}.*rgb0(k0,1);
      im2R0(:,:,2) = imPattern{indImPattern}.*rgb0(k0,2);
      im2R0(:,:,3) = imPattern{indImPattern}.*rgb0(k0,3);
      im2R1(:,:,1) = imPattern{indImPattern}.*rgb1(k0,1);
      im2R1(:,:,2) = imPattern{indImPattern}.*rgb1(k0,2);
      im2R1(:,:,3) = imPattern{indImPattern}.*rgb1(k0,3);
      %wn=cwin0(img1, 'Stereo', cf, rc00, window1, window2);
      [iLf0 iRf0]=cwin3(im2R0, im2R0, cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(k0));
      % zaber(name_map('rotation')).move_deg(dgs(k0)); %%-6400

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
                  if keyCode(KbName('RightArrow')) | keyCode(KbName('5'))
                      opt_chk = 1;
                      %end
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
      
      if scene.enable_tcp; send_tcp0fiatAcu(tcp_socket, 1, k0, vs); end; t0(k0,:)=clock;
      
      snd(1000, 0.2); pause(0.5);
      for i = 1:floor(length(sinValues)/2)
         opto(name_map('l_disp')).control.setFocalPower(power_dispL-sinValues(i));
         opto(name_map('r_disp')).control.setFocalPower(power_dispR-sinValues(i));
         pause(tIntervalStm);
      end
      tChange1(k0,:) = clock;
      [iLf1 iRf1]=cwin3(im2R1, im2R1, cf, rc00, window1, window2);
      tChange2(k0,:) = clock;
      for i = (floor(length(sinValues)/2)+1):length(sinValues)
         opto(name_map('l_disp')).control.setFocalPower(power_dispL-sinValues(i));
         opto(name_map('r_disp')).control.setFocalPower(power_dispR-sinValues(i));
         pause(tIntervalStm);
      end      
      snd(1000, 0.1); pause(0.2);
      snd(1000, 0.1); pause(0.1);
     
      if scene.enable_tcp; send_tcp0fiatAcu(tcp_socket, 0, k0, vs); end; %stage) 0stop 1record
      tRealEnd(k0,:) = clock;
      %pause(3);
      %wn=cwin0(img0, 'Stereo', cf, rc00, window1, window2);
      [iLf0 iRf0]=cwin3(zeros(size(im2R0)), zeros(size(im2R0)), cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(k0));
    %           zaber(name_map('rotation')).move_deg(-3); %%-6400
    %  zaber(name_map('rotation')).move_deg(dgs0); %%-6400

      t1(k0,:)=clock;

      %snd(2000, 0.25);  pause(2.75);
      pause(1);

      %if scene.enable_tcp; send_tcp0(scene, 0); end %stage) 0stop 1record
      t2(k0,:)=clock;

    scene.trial_num=k0;
end

% if scene.enable_tcp
%     clear tcp_socket;
% end
AFCp.v1=power_dispR
AFCp.t3=cat(3, t0, t1, t2,tChange1,tChange2,tRealEnd);
AFCp.dgs=dgs;
AFCp.sinValues = sinValuesAll;
AFCp.indImPatternAll = indImPatternAll;
