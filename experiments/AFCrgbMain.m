function AFCp=AFCrgbMain(imPattern, rgb0, meanv0, sr, window1, window2, vs)

% meanv0 is vector of optical distances
%        control: [1×1 Zaber.AsciiDevice]
%     unit_scale: 2.1333e+03
%           move: @(x)zaber(name_map(ident)).control.moveabsolute(x)
%       move_deg: @(x)zaber(name_map(ident)).control.moveabsolute(x*zaber(name_map(ident)).unit_scale)
%        move_mm: []

global sz cf rc00 name_map zaber opto log

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

t0=zeros(length(meanv0), 6); t1=t0; t2=t0; tRealEnd = t0;
disp('ready to start');  KbWait([], 2); 
indImPatternAll = [];
for k0=1:size(meanv0,1)
      if length(imPattern)>1
        indImPattern = randsample(1:length(imPattern),1);
      else
        indImPattern = 1;
      end
      indImPatternAll(end+1,:) = indImPattern;
      im2R0 = [];
      im2R0(:,:,1) = imPattern{indImPattern}.*rgb0(k0,1);
      im2R0(:,:,2) = imPattern{indImPattern}.*rgb0(k0,2);
      im2R0(:,:,3) = imPattern{indImPattern}.*rgb0(k0,3);
      %wn=cwin0(img1, 'Stereo', cf, rc00, window1, window2);
      [iLf0 iRf0]=cwin3(im2R0, im2R0, cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(k0));

      %disp( n2s(v0(k0)));        
      fprintf('TRL= %f, L = %f  , R = %f , DEG = %f, Demand = %f\n' ,k0, opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power, (zaber(name_map('rotation')).control.getposition)./2.1333E3, meanv0(k0) );

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
                  % if keyCode(KbName('RightArrow')) | keyCode(12)
                  if keyCode(KbName('RightArrow')) | keyCode(KbName('5')) | keyCode(12)
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
      snd(1000, 0.2); 
      pause(3.0);
      if scene.enable_tcp; send_tcp0fiatAcu(tcp_socket, 0, k0, vs); end; %stage) 0stop 1record
      tRealEnd(k0,:) = clock;

      [iLf0 iRf0]=cwin3(zeros(size(im2R0)), zeros(size(im2R0)), cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanv0(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanv0(k0));
    %           zaber(name_map('rotation')).move_deg(-3); %%-6400

      t1(k0,:)=clock;
      
      if k0<size(meanv0,1) && ~isequal(rgb0(k0+1,:),rgb0(k0,:))
         pause(3);
      else
         pause(0.5);
      end

      %if scene.enable_tcp; send_tcp0(scene, 0); end %stage) 0stop 1record
      t2(k0,:)=clock;

    scene.trial_num=k0;
end

% if scene.enable_tcp
%     clear tcp_socket;
% end
AFCp.v1=power_dispR
AFCp.t3=cat(3, t0, t1, t2,tRealEnd);
AFCp.indImPatternAll = indImPatternAll;
