%%211102 DSP_IMG

ey = 'R';
[window1, window2, vbl0]=strt_psych0(screenNumber-2, screenNumber-1, 0);
power_dispL_min = 7;
power_dispL_max = 16;
power_dispR_min = 7;
power_dispR_max = 16.4;
adjustIncrement = 0.1;
rgbTest = [1 1 1];
rgbTest = round(rgbTest.*255);

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
opto(name_map('r_disp')).control.setFocalPower(14.4+sr(2));% -dmnd(k0));
opto(name_map('r_disp')).control.getFocalPower.focal_power

opto(name_map('l_disp')).control.setFocalPower(14+sr(1));%-dmnd(k0));
opto(name_map('l_disp')).control.getFocalPower.focal_power

bTexture = false;
testim = [];
testim(:,:,1) = ones([900 900]).*rgbTest(1);
testim(:,:,2) = ones([900 900]).*rgbTest(2);
testim(:,:,3) = ones([900 900]).*rgbTest(3);

[iLf iRf]=cwin3(imread("black.png"), testim , cf, rc00, window2, window1);

power_dispL = 14;
power_dispR = 13.9;

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
            elseif (keyCode(KbName('DownArrow')) | keyCode(KbName('5')))
                [iLf iRf]=cwin3(imread("black.png"), testim , cf, rc00, window2, window1);  
                display('Refresh texture');
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
            
            fprintf('Display power: L = %f  , R = %f , Optical Distance R = %f D \n',power_dispL, power_dispR, 0.87.*(13.935-power_dispR));

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

sca   
