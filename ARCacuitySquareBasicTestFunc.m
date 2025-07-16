function AFCp = ARCacuitySquareBasicTestFunc(imPattern,rgb,meanFocstmOptDst,frqCpd, contrast, window1, window2, trlPerLvl)

global cf rc00 name_map zaber opto log

rgbAll = [];
meanFocstmOptDstAll = [];
contrastAll = [];
maskBrightness = 0;
gammaR = 2.5;
gammaG = 2.7;
gammaB = 2.3;
bAccStimEqualsAcuStim = false;
powerDiscrepancy = 0;

for i = 1:size(rgb,1)
   for j = 1:size(contrast,2)
       for l = 1:length(meanFocstmOptDst)
           for k = 1:trlPerLvl
               rgbAll(end+1,:) = rgb(i,:);
               contrastAll(end+1,:) = contrast(j);
               meanFocstmOptDstAll(end+1,:) = meanFocstmOptDst(l);
           end
       end
   end
end
indScramble = randperm(length(contrastAll))';
rgbAll = rgbAll(indScramble,:);
contrastAll = contrastAll(indScramble);
meanFocstmOptDstAll = meanFocstmOptDstAll(indScramble);
rgbAll(end+1,:) = [0 0 0];
contrastAll(end+1,:) = 1;
meanFocstmOptDstAll(end+1,:) = 3;
meanFocstmOptDstAllOrig = meanFocstmOptDstAll;

% 1 = 0°, 2 = 90°, 3 = 180°, 4 = 270° 
stimOrientation = ceil(rand(size(contrastAll))*2);

power_dispR=14.3; %starting display power
power_dispL=14; %starting display power
opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(1));
opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(1));

if length(imPattern)>1
   indImPattern = randsample(1:length(imPattern),1);
else
   indImPattern = 1;
end

im2R0 = [];
if bAccStimEqualsAcuStim
    accStim = ARC2Dgabor(smpPos(512,512),[],0,0,[8 24 40 56],[0.9 0.9/3 0.9/5 0.9/7],0,90,0.3,0.3,[rgb(1,1)^gammaR rgb(1,2)^gammaG rgb(1,3)^gammaB],1,1,0,0);
    accStim(:,:,1) = accStim(:,:,1).^(1/gammaR);
    accStim(:,:,2) = accStim(:,:,2).^(1/gammaG);
    accStim(:,:,3) = accStim(:,:,3).^(1/gammaB);    
    im2R0(:,:,1) = accStim(:,:,1).*255;
    im2R0(:,:,2) = accStim(:,:,2).*255;
    im2R0(:,:,3) = accStim(:,:,3).*255;
else
   im2R0(:,:,1) = imPattern{indImPattern}.*rgb(1,1);
   im2R0(:,:,2) = imPattern{indImPattern}.*rgb(1,2);
   im2R0(:,:,3) = imPattern{indImPattern}.*rgb(1,3);
end

cwin3(im2R0, im2R0, cf, rc00, window1, window2);

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

t0=zeros(length(contrastAll), 6); t1=t0; t2=t0; tChange1 = t0; tChange2 = t0; tRealEnd = t0;
% stage) 0stop 1record figure this out with Steve
disp('ready to start');  KbWait([], 2); 

rspAcu = [];
for k0=1:length(contrastAll)
      if length(imPattern)>1
        indImPattern = randsample(1:length(imPattern),1);
      else
        indImPattern = 1;
      end
      imPatternTmp = imPattern{indImPattern};
      if ~bAccStimEqualsAcuStim
          im2R0 = [];
          im2R0(:,:,1) = imPatternTmp.*rgbAll(k0,1);
          im2R0(:,:,2) = imPatternTmp.*rgbAll(k0,2);
          im2R0(:,:,3) = imPatternTmp.*rgbAll(k0,3);
      end
      blackStim = zeros(size(im2R0));
      acuStimOrig1 = ARC2Dgabor(smpPos(256,256),[],0,0,[frqCpd 3*frqCpd 5*frqCpd 7*frqCpd],[contrastAll(k0) contrastAll(k0)/3 contrastAll(k0)/5 contrastAll(k0)/7],-15,90,0.2,0.2,[rgbAll(k0,1)^gammaR rgbAll(k0,2)^gammaG rgbAll(k0,3)^gammaB],1,1,0,0);
      acuStimOrig1(:,:,1) = acuStimOrig1(:,:,1).^(1/gammaR);
      acuStimOrig1(:,:,2) = acuStimOrig1(:,:,2).^(1/gammaG);
      acuStimOrig1(:,:,3) = acuStimOrig1(:,:,3).^(1/gammaB);
      acuStimOrig2 = ARC2Dgabor(smpPos(256,256),[],0,0,[frqCpd 3*frqCpd 5*frqCpd 7*frqCpd],[contrastAll(k0) contrastAll(k0)/3 contrastAll(k0)/5 contrastAll(k0)/7],15,90,0.2,0.2,[rgbAll(k0,1)^gammaR rgbAll(k0,2)^gammaG rgbAll(k0,3)^gammaB],1,1,0,0);
      acuStimOrig2(:,:,1) = acuStimOrig2(:,:,1).^(1/gammaR);
      acuStimOrig2(:,:,2) = acuStimOrig2(:,:,2).^(1/gammaG);
      acuStimOrig2(:,:,3) = acuStimOrig2(:,:,3).^(1/gammaB);      
      acuStimOrig(:,:,:,1) = acuStimOrig1;
      acuStimOrig(:,:,:,2) = acuStimOrig2;
      acuStim = acuStimOrig.*255;    
      noiseStim = maskBrightness.*repmat(rand([5 5 1]),[1 1 3]);
      noiseStim = imresize(noiseStim,[100 100],'nearest');      
      cwin3(im2R0, im2R0, cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0));

      fprintf('TRL= %f, L = %f  , R = %f , DEG = %f, Demand = %f\n' ,k0, opto(name_map('l_disp')).control.getFocalPower.focal_power, opto(name_map('r_disp')).control.getFocalPower.focal_power, (zaber(name_map('rotation')).control.getposition)./2.1333E3, meanFocstmOptDstAll(k0) );

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
                  if keyCode(KbName('RightArrow')) | keyCode(KbName('6'))
                      opt_chk = 1;
                      if k0>1 rspAcu(k0-1) = 1; end
                  elseif keyCode(KbName('LeftArrow')) | keyCode(KbName('4'))
                      opt_chk = 1;
                      if k0>1 rspAcu(k0-1) = 2; end                 
                  elseif keyCode(KbName('UpArrow')) | keyCode(KbName('8'))
                      powerDiscrepancy = powerDiscrepancy+0.25/0.816;
                      disp(['Power discrepancy = ' num2str(powerDiscrepancy)]);
                      meanFocstmOptDstAll = meanFocstmOptDstAllOrig+powerDiscrepancy;
                      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0));
                      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0));                      
                  elseif keyCode(KbName('DownArrow')) | keyCode(KbName('2'))
                      powerDiscrepancy = powerDiscrepancy-0.25/0.816;            
                      disp(['Power discrepancy = ' num2str(powerDiscrepancy)]);
                      meanFocstmOptDstAll = meanFocstmOptDstAllOrig+powerDiscrepancy;
                      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0));
                      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0));                            
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
      if k0==length(contrastAll)
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
      % opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0));
      % opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0));
      cwin3(im2R0, im2R0, cf, rc00, window1, window2);
      pause(1);
      cwin3(blackStim, blackStim, cf, rc00, window1, window2);
      tChange1(k0,:) = clock;
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0));      
      if scene.enable_tcp; send_tcp0(scene, 1); end; t1(k0,:)=clock;
      pause(0.1);      
      cwin3(squeeze(acuStim(:,:,:,stimOrientation(k0))), squeeze(acuStim(:,:,:,stimOrientation(k0))), cf, rc00, window1, window2);
      tChange2(k0,:) = clock;
      pause(0.100);
      cwin3(noiseStim, noiseStim, cf, rc00, window1, window2);
      opto(name_map('l_disp')).control.setFocalPower(power_dispL-meanFocstmOptDstAll(k0));
      opto(name_map('r_disp')).control.setFocalPower(power_dispR-meanFocstmOptDstAll(k0));      
      pause(0.15);
      cwin3(blackStim, blackStim, cf, rc00, window1, window2);
      if scene.enable_tcp; send_tcp0(scene, 0); end %stage) 0stop 1record
      t2(k0,:)=clock;
      % opto(name_map('l_disp')).control.setFocalPower(power_dispL);
      % opto(name_map('r_disp')).control.setFocalPower(power_dispR);

      pause(0.2);
      scene.trial_num=k0;
      tRealEnd(k0,:) = clock;
      display(['k0 = ' num2str(k0) ', rsp length = ' num2str(length(rspAcu))]);
end

if scene.enable_tcp; fclose(scene.tcp_socket); end
AFCp.v1=power_dispR;
t3=cat(3, t0, t1, t2,tChange1,tChange2,tRealEnd);
AFCp.t3 = t3(1:end-1,:,:);
AFCp.rgb = rgbAll(1:end-1,:);
AFCp.meanFocstmOptDst = meanFocstmOptDstAll(1:end-1);
AFCp.frqCpd = frqCpd;
AFCp.contrast = contrastAll(1:end-1);
AFCp.rspAcu = rspAcu;
AFCp.stimOrientation = stimOrientation(1:end-1);
AFCp.im2R0 = im2R0;
AFCp.acuStim = acuStim;
