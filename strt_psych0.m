function [window1, window2 vbl]=strt_psych0(screenNumber1, screenNumber2, c1);
global ek uk dk lk rk sp kntr n0m n1m n2m n0d n1d n2d black white grey windowRect ifi  screenXpixels screenYpixels xCenter yCenter

psych0(2)

black1 = BlackIndex(screenNumber1);
white1 = WhiteIndex(screenNumber1);
grey1 = white1 / 2;

black2 = BlackIndex(screenNumber2);
white2 = WhiteIndex(screenNumber2);
grey2 = white2 / 2;

black=black1;
white=white1;
grey=grey1;

if c1==0;
    clr0=black;
elseif c1==1;
    clr0=white;
elseif c1==2;
    clr0=grey;
end

% switch ex
%     case 'Tumbling E'
%         clr0=white;
%     otherwise
%         clr0=black
% end
    


% Open an on screen window and color it grey
[window1, windowRect1] = PsychImaging('OpenWindow', screenNumber1, clr0);

[window2, windowRect2] = PsychImaging('OpenWindow', screenNumber2, clr0);

windowRect=windowRect1;

% Set the blend funciton for the screen
Screen('BlendFunction', window1, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Set the blend funciton for the screen
Screen('BlendFunction', window2, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');


% Query the frame duration
ifi1 = Screen('GetFlipInterval', window1);

ifi2 = Screen('GetFlipInterval', window2);

ifi=ifi1;

% Get the size of the on screen window in pixels
% For help see: Screen WindowSize?
[screenXpixels1, screenYpixels1] = Screen('WindowSize', window1);

[screenXpixels2, screenYpixels2] = Screen('WindowSize', window2);

screenXpixels=screenXpixels1;
screenYpixels=screenYpixels1;

% Get the centre coordinate of the window in pixels
% For help see: help RectCenter
[xCenter1, yCenter1] = RectCenter(windowRect1);

[xCenter2, yCenter2] = RectCenter(windowRect2);

xCenter=xCenter1;
yCenter=yCenter1;

% The avaliable keys to press
ek = KbName('ESCAPE');
uk = KbName('UpArrow');
dk = KbName('DownArrow');
lk = KbName('LeftArrow');
rk = KbName('RightArrow');
sp = KbName('space'); %'Return' if not good enough
n0m = KbName('0)'); %'
n1m = KbName('1!'); %'
n2m = KbName('2@'); %'
n0d = KbName('0');
n1d = KbName('1'); %'
n2d = KbName('2'); %'
kntr = KbName('Return'); %'

%KbName('KeyNames')
vbl(1) = Screen('Flip', window1);
vbl(2) = Screen('Flip', window2);

% function strt_psych(screenNumber)
% global ek uk dk lk rk st black white grey screenXpixels screenYpixels xCenter yCenter window ifi
% 
% psych0(2)
% 
% black = BlackIndex(screenNumber);
% white = WhiteIndex(screenNumber);
% grey = white / 2;
% 
% % Open an on screen window and color it grey
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
% 
% % Set the blend funciton for the screen
% Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
% 
% % Query the frame duration
% ifi = Screen('GetFlipInterval', window);
% 
% % Get the size of the on screen window in pixels
% % For help see: Screen WindowSize?
% [screenXpixels, screenYpixels] = Screen('WindowSize', window);
% 
% % Get the centre coordinate of the window in pixels
% % For help see: help RectCenter
% [xCenter, yCenter] = RectCenter(windowRect);
% 
% % The avaliable keys to press
% ek = KbName('ESCAPE');
% uk = KbName('UpArrow');
% dk = KbName('DownArrow');
% lk = KbName('LeftArrow');
% rk = KbName('RightArrow');
