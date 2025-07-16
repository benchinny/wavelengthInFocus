%%211102 DSP_IMG

ey = 'R';

power_dispL_min = 7;
power_dispL_max = 16;
power_dispR_min = 7;
power_dispR_max = 16.4;
adjustIncrement = 0.12255;
bRecord = 1;
clrIndAll = repmat([1 2 3]',[40 1]);
clrIndAll = clrIndAll(randperm(length(clrIndAll)));
clrIndCounter = 1;

if bRecord && ~ismember('tcp_socket', who('global'))
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
end

[window1, window2, vbl0]=strt_psych0(screenNumber-2, screenNumber-1, 0);

%%input a output b
cf=ones(3,2);

%initial  power=13.5 AR    
%LCAim='texture0_1080_newfill_malt.png';    
%im2R='black.png';       
% im2L='texture0_1080_newfill_malt.png';
% im2L = 'G:\My Drive\exp_bvams\code_repo\imgs\TCA_r540_k120_b108_cw5_sbm7.png';  
%im2L='vrn10_G_sd1.png';

deg=-3;            
zaber(name_map('rotation')).move_deg(deg); %%-6400            
zaber(name_map('rotation')).control.getposition  

if exist('sr') ~=  1; sr=[0 0]; end
%         dmnd=[-0.5:0.5:3]
opto(name_map('r_disp')).control.setFocalPower(15.3+sr(2));% -dmnd(k0));
opto(name_map('r_disp')).control.getFocalPower.focal_power

opto(name_map('l_disp')).control.setFocalPower(14+sr(1));%-dmnd(k0));
opto(name_map('l_disp')).control.getFocalPower.focal_power

bTexture = false;
if bTexture
    % ----LOADING IMAGE----------
    % im2L='texture0_nrm_rgb.png';
    im2L='texture0_1080_newfill_malt.png';
    % im2L = 'TCA_r540_k120_b120_cw5.png';
    % im2L = 'testEresized.png';
    im2R=im2L;
    testim = imread(im2R);
%    testim = 255.*ones(size(testim));
    % ---------------------------
    % ----MAKING GABOR-----------
    % testim = ARC2Dgabor(smpPos(512,512),[],0,0,30,1,-45,0,0.16,0.5,[0.25 0 0],1,1,0,0);
    % testim(:,:,1) = testim(:,:,1).^(1/2.4);
    % testim(:,:,3) = testim(:,:,3).^(1/2.2);    
    % testim = testim.*255;
    % ---------------------------
    % ----MAKING TUMBLING E------
%     acuStimOrig = imread('testEresized.png');
%     indGdAcu = acuStimOrig>0;
%     acuStimOrig(indGdAcu) = 255;
%     acuStimOrig(~indGdAcu) = 0;
%     acuStim = [];
%     acuStimTmp = imresize(acuStimOrig,200.*[1 1]);
%     acuStimTmpRGB = [];
%     acuStimTmpRGB(:,:,1) = acuStimTmp.*0.56;
%     acuStimTmpRGB(:,:,2) = acuStimTmp.*0;
%     acuStimTmpRGB(:,:,3) = acuStimTmp.*0;   
%     testim = acuStimTmpRGB;
    % ---------------------------
else
    clrInd = clrIndAll(1);
    im1 = imread('H:\Shared drives\CIVO_BVAMS\stimuli\word_image_01.png');
    im1(im1>0) = 255;
    im1 = flipud(im1);   
    imPatternTmp = squeeze(im1(:,:,3));
    imPatternTmp = [zeros([30 size(imPatternTmp,2)]); imPatternTmp; zeros([30 size(imPatternTmp,2)])];
    imPatternTmp = [zeros([size(imPatternTmp,1) 30]) imPatternTmp zeros([size(imPatternTmp,1) 30])];  

    imPatternColor = [];
    imPatternColor(:,:,1) = round(imPatternTmp.*0.569);
    imPatternColor(:,:,2) = round(imPatternTmp.*0.432);
    imPatternColor(:,:,3) = imPatternTmp;

    testim = zeros([size(imPatternTmp) 3]);
    testim(:,:,clrInd) = squeeze(imPatternColor(:,:,clrInd));
    display(['Color index is ' num2str(clrInd)]);
end
[iLf iRf]=cwin3(imread("black.png"), testim , cf, rc00, window2, window1);

power_dispL = 14;
power_dispR = 15.3;
powerDispRall = [];

rightTrombonePowerNear = opto(name_map('r_t_near')).control.getFocalPower.focal_power;
rightTrombonePowerFar = opto(name_map('r_t_far')).control.getFocalPower.focal_power;
rightTrombonePosition = zaber(name_map('r_trombone')).control.getposition;

fprintf('Trombone Near R = %f  , Trombone Far R = %f\n',rightTrombonePowerNear, rightTrombonePowerFar);

fprintf('\n');

addpath(genpath(fullfile('toolboxes')));
KbName('UnifyKeyNames');
KbWait([], 2); 

%% Control loop
ListenChar(2);
try
    opt_chk=0;
    while opt_chk==0
        [ keyIsDown, keyTime, keyCode ] = KbCheck;
        if keyIsDown
            if keyCode(KbName('RightArrow')) | keyCode(KbName('6'))
                if ey(1)=='R'
                    power_dispR=power_dispR+adjustIncrement;
                elseif ey(1)=='L'
                    power_dispL=power_dispL+adjustIncrement;
                else 
                    power_dispR=power_dispR+adjustIncrement;
                    power_dispL=power_dispL+adjustIncrement;
                end
                %end
            elseif keyCode(KbName('LeftArrow')) | keyCode(KbName('4'))
                if ey(1)=='R'
                    power_dispR=power_dispR-adjustIncrement;
                elseif ey(1)=='L'
                    power_dispL=power_dispL-adjustIncrement;
                else
                    power_dispR=power_dispR-adjustIncrement;
                    power_dispL=power_dispL-adjustIncrement;
                end
            elseif ~bTexture & (keyCode(KbName('DownArrow')) | keyCode(KbName('5')) | keyCode(12))
                if bRecord
                    % RECORD CURRENT WAVEFRONT
                    send_tcp0fiatAcu(tcp_socket, 1, clrInd, vs);
                    pause(0.25);
                    send_tcp0fiatAcu(tcp_socket, 0, clrInd, vs);
                end
                
                % SAVE CURRENT OPTOTUNE POWER AND COLOR INDEX
                powerDispRall(end+1) = power_dispR;
                
                % ZERO OUT STIMULUS AT CURRENT COLOR INDEX
                testim(:,:,clrInd) = zeros([size(testim,1) size(testim,2)]);
                
                clrIndCounter = clrIndCounter+1;

                clrInd = clrIndAll(clrIndCounter); % NEW COLOR INDEX
                testim(:,:,clrInd) = squeeze(imPatternColor(:,:,clrInd));
                display(['Color index is ' num2str(clrInd)]);
                power_dispR = 15.3;
                [iLf iRf]=cwin3(imread("black.png"), testim , cf, rc00, window2, window1);
                if clrIndCounter>60
                    opt_chk = 1;
                end
            elseif keyCode(KbName('Return')) %| keyCode(KbName('Return'))
                opt_chk=1;    
            else
                disp('WRONG KEY'); snd(100, 0.25);
            end
            
            if power_dispL < power_dispL_min
                power_dispL = power_dispL_min;
            elseif power_dispL > power_dispL_max;
                power_dispL = power_dispL_max;
            end

            if power_dispR < power_dispR_min
                power_dispR = power_dispR_min;
            elseif power_dispR > power_dispR_max;
                power_dispR = power_dispR_max;
            end            

            opto(name_map('l_disp')).control.setFocalPower(power_dispL);
            opto(name_map('r_disp')).control.setFocalPower(power_dispR);
            
            fprintf('Display power: L = %f  , R = %f , Optical Distance R = %f D, Trial %f \n',power_dispL, power_dispR, 0.816.*(15.3-power_dispR),clrIndCounter);

            fprintf('\n');
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

KbWait([], 2);
[iLf iRf]=cwin3(imread("black.png"), imread("black.png") , cf, rc00, window2, window1);
clrIndAll = clrIndAll(1:clrIndCounter);
filePath = 'H:\Shared drives\CIVO_BVAMS\data\ARC\';
save([filePath 'LCAfileS' num2str(sn-1000) '.mat'],'powerDispRall','clrIndAll');
sca   
% if bRecord
%     clear tcp_socket;
% end
