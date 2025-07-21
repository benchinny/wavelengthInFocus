function AFCp = ARCacuitySquareFunc(imPattern,rgb,meanFocstmOptDst,focStmOptDstIncr, window1, window2, trlPerLvl,vs,frqCpd, contrast)

global cf rc00 name_map zaber opto log

rgbAll = [];
meanFocstmOptDstAll = [];
focStmOptDstIncrAll = [];
indScramble = [];
maskBrightness = 0;
maskSize = [100 100];
gammaR = 2.5;
gammaG = 2.7;
gammaB = 2.3;
nRec = 1;

% CAREFUL ATTEMPT TO BLOCK CONDITIONS SO EACH OPTICAL DISTANCE INCREMENT IS
% PRESENTED ONCE PER BLOCK
for i = 1:size(rgb,1)
   for j = 1:length(meanFocstmOptDst)
       for k = 1:length(focStmOptDstIncr)
           rgbAll(end+1,:) = rgb(i,:);
           meanFocstmOptDstAll(end+1,:) = meanFocstmOptDst(j);
           focStmOptDstIncrAll(end+1,:) = focStmOptDstIncr(k);
       end
       for l = 1:trlPerLvl
          indScramble = [indScramble; randperm(length(focStmOptDstIncr))'];
       end
   end
end
% RANDOMIZING TRIALS
indScramble = indScramble+imresize(length(focStmOptDstIncr).*[0:(trlPerLvl*size(rgb,1)*length(meanFocstmOptDst)-1)]',size(indScramble),'nearest');
rgbAll = repmat(rgbAll,[trlPerLvl 1]);
meanFocstmOptDstAll = repmat(meanFocstmOptDstAll,[trlPerLvl 1]);
focStmOptDstIncrAll = repmat(focStmOptDstIncrAll,[trlPerLvl 1]);
stimSizePixAll = 10.*ones(size(focStmOptDstIncrAll));
offsetXall = 5.*ones(size(focStmOptDstIncrAll));
offsetYall = 10.*ones(size(focStmOptDstIncrAll));
sizeTotal = 100;
stimCtr = 50;
rgbAll = rgbAll(indScramble,:);
meanFocstmOptDstAll = meanFocstmOptDstAll(indScramble);
focStmOptDstIncrAll = focStmOptDstIncrAll(indScramble);
stimSizePixAll = stimSizePixAll(indScramble);

% ADD DUMMY TRIAL RIGHT AT THE END (PECULIAR TO WAY CODE IS WRITTEN)
rgbAll(end+1,:) = [0 0 0];
focStmOptDstIncrAll(end+1,:) = 0;
meanFocstmOptDstAll(end+1,:) = 3;
stimSizePixAll(end+1,:) = 10;

% 1 = 0°, 2 = 90°, 3 = 180°, 4 = 270° 
stimOrientation = ceil(rand(size(focStmOptDstIncrAll))*2);

power_dispR=14.3; %starting display power
power_dispL=14; %starting display power
opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(1));
opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(1));

if length(imPattern)>1
   indImPattern = randsample(1:length(imPattern),1);
else
   indImPattern = 1;
end
im2R0(:,:,1) = imPattern{indImPattern}.*rgb(1,1);
im2R0(:,:,2) = imPattern{indImPattern}.*rgb(1,2);
im2R0(:,:,3) = imPattern{indImPattern}.*rgb(1,3);

cwin3(im2R0, im2R0, cf, rc00, window1, window2);

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

t0=zeros(length(focStmOptDstIncrAll), 6); t1=t0; t2=t0; tChange1 = t0; tChange2 = t0; tRealEnd = t0;
% stage) 0stop 1record figure this out with Steve
disp('ready to start');  KbWait([], 2); 

rspAcu = [];
indImPatternAll = [];
for k0=1:length(focStmOptDstIncrAll)
      if length(imPattern)>1
        indImPattern = randsample(1:length(imPattern),1);
      else
        indImPattern = 1;
      end
      indImPatternAll(end+1,:) = indImPattern; 
      im2R0 = [];
      im2R0(:,:,1) = imPattern{indImPattern}.*rgbAll(k0,1);
      im2R0(:,:,2) = imPattern{indImPattern}.*rgbAll(k0,2);
      im2R0(:,:,3) = imPattern{indImPattern}.*rgbAll(k0,3);
      blackStim = zeros(size(im2R0));
      acuStimOrig1 = ARC2Dgabor(smpPos(256,256),[],0,0,[frqCpd 3*frqCpd 5*frqCpd 7*frqCpd],[contrast contrast/3 contrast/5 contrast/7],-15,90,0.2,0.2,[rgbAll(k0,1)^gammaR rgbAll(k0,2)^gammaG rgbAll(k0,3)^gammaB],1,1,0,0);
      acuStimOrig1(:,:,1) = acuStimOrig1(:,:,1).^(1/gammaR);
      acuStimOrig1(:,:,3) = acuStimOrig1(:,:,3).^(1/gammaB);
      acuStimOrig2 = ARC2Dgabor(smpPos(256,256),[],0,0,[frqCpd 3*frqCpd 5*frqCpd 7*frqCpd],[contrast contrast/3 contrast/5 contrast/7],15,90,0.2,0.2,[rgbAll(k0,1)^gammaR rgbAll(k0,2)^gammaG rgbAll(k0,3)^gammaB],1,1,0,0);
      acuStimOrig2(:,:,1) = acuStimOrig2(:,:,1).^(1/gammaR);
      acuStimOrig2(:,:,3) = acuStimOrig2(:,:,3).^(1/gammaB);      
      acuStimOrig(:,:,:,1) = acuStimOrig1;
      acuStimOrig(:,:,:,2) = acuStimOrig2;
      acuStim = acuStimOrig.*255;
      noiseStim = maskBrightness.*repmat(rand([5 5 1]),[1 1 3]);
      noiseStim = imresize(noiseStim,maskSize,'nearest');
      cwin3(im2R0, im2R0, cf, rc00, window1, window2);
      if mod(k0,length(focStmOptDstIncr))==1
         snd(250, 0.25);
         snd(500, 0.25);
      else
         snd(250, 0.25); %pause(2.75);
      end      

      fprintf('TRL= %f, L = %f  , R = %f , DEG = %f, Demand = %f\n' ,k0, opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power, (zaber(name_map('rotation')).control.getposition)./2.1333E3, meanFocstmOptDstAll(k0) );
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
                  if keyCode(KbName('RightArrow')) | keyCode(KbName('6'))
                      opt_chk = 1;
                      if k0>1 rspAcu(k0-1) = 1; end
                  elseif keyCode(KbName('LeftArrow')) | keyCode(KbName('4'))
                      opt_chk = 1;
                      if k0>1 rspAcu(k0-1) = 2; end
                  elseif keyCode(KbName('Escape')) %| keyCode(KbName('Return'))
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
      if k0==length(focStmOptDstIncrAll)
          break;
      end      
      t0(k0,:)=clock;
      if k0>1 && rspAcu(k0-1)==stimOrientation(k0-1)
         snd(1000, 0.2);
      elseif k0>1 && rspAcu(k0-1)~=stimOrientation(k0-1)
         snd(200, 0.2);
      elseif k0==1
         snd(1000,0.2);
      end
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0));      
      % opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0));
      % opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0));
      cwin3(im2R0, im2R0, cf, rc00, window1, window2);
      pause(2);
      t1(k0,:)=clock;
      cwin3(blackStim, blackStim, cf, rc00, window1, window2);
      tChange1(k0,:) = clock;
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0)-focStmOptDstIncrAll(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0)-focStmOptDstIncrAll(k0));      
      if scene.enable_tcp
          send_tcp0fiatAcu(tcp_socket, 1, k0, vs); 
      end       
      pause(0.15);
      cwin3(squeeze(acuStim(:,:,:,stimOrientation(k0))), squeeze(acuStim(:,:,:,stimOrientation(k0))), cf, rc00, window1, window2);
      tChange2(k0,:) = clock;
      pause(0.10);
      cwin3(noiseStim, noiseStim, cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0));      
      pause(0.15);
      if scene.enable_tcp
          send_tcp0fiatAcu(tcp_socket, 0, k0, vs)
      end %stage) 0stop 1record
      cwin3(blackStim, blackStim, cf, rc00, window1, window2);
      t2(k0,:)=clock;

      pause(0.2);
      scene.trial_num=k0;
      tRealEnd(k0,:) = clock;
end

% if scene.enable_tcp
%     clear tcp_socket;
% end
AFCp.v1=power_dispR;
t3=cat(3, t0, t1, t2,tChange1,tChange2,tRealEnd);
AFCp.t3 = t3(1:end-1,:,:);
AFCp.rgb = rgbAll(1:end-1,:);
AFCp.meanFocstmOptDst = meanFocstmOptDstAll(1:end-1);
AFCp.focStmOptDstIncr = focStmOptDstIncrAll(1:end-1);
AFCp.rspAcu = rspAcu;
AFCp.stimOrientation = stimOrientation(1:end-1);
AFCp.im2R0 = im2R0;
AFCp.stimSizePixAll = stimSizePixAll(1:end-1);
AFCp.contrast = contrast;
AFCp.frqCpd = frqCpd;
AFCp.indImPatternAll = indImPatternAll;
