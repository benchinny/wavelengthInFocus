%% Version 210510
%clr
global ek uk dk lk rk st snt  black white grey screenXpixels screenYpixels xCenter yCenter ifi zbr cf Ec sz E10c rc00 name_map l_trombone_f r_trombone_f l_opto_f r_opto_f  enable_optotunes enable_trombones zaber opto ey cntrst log
% edit JJVC_stro
% addpath([pwd '\JnJ']); addpath([pwd '\fcns']); addpath([pwd '\imgs']); addpath([pwd '\apps']); addpath([pwd '\toolboxes']); addpath([pwd '\data']); addpath([pwd '\psignifit-master']); addpath([pwd '\data\System TCA\Data Processing']);  addpath(genpath(fullfile('toolboxes')));
% addpath([pwd '\JnJ']); addpath([pwd '\fcns']); addpath([pwd '\imgs']); addpath([pwd '\apps']); addpath([pwd '\toolboxes']); addpath([pwd '\psignifit-master']); addpath([pwd '\data\System TCA\Data Processing']);  addpath(genpath(fullfile('toolboxes')));
addpath([pwd '\fcns']); addpath('H:\Shared drives\CIVO_BVAMS\code\imgs\'); addpath([pwd '\apps']); addpath([pwd '\toolboxes']); addpath([pwd '\psignifit-master']); addpath(['H:\Shared drives\CIVO_BVAMS\data\System TCA\Data Processing\']);  addpath(genpath(fullfile('toolboxes'))); ARCtools;
cls
% 0.00384deg/pixels,  260.417 pixels/degree, 4.34 pixels/arcmin,  0.7234 pixels/arcsec
c0=0.00384; %deg/pixels
%display =10 degrees
%FOV= 5 degrees

%BVAMS_gui; %run input gui
% ey='Right'; %eye Right Left Binocular
ey='Binc'; %eye Right Left Binocular
% magn=0.8;
magn=1;
ACL=0; %ACL attached=1 otherwise=0  %!!!NOTE!!! NEED TO SET VALUE OF ACL correctly for findz0 to work correctly!
sn=1001; %subject  number 10004MB
vs=2; % visit number
load('JnJrand211130', 'ETMm', 'AFCm', 'TCAm')
%sr=[0 0]; %[L R] spherical refraction
%subject number 10003AR
%subject number 10007 Nadav
%subject number 10005 Jim
%subject number 10009 Fabio
%sn=77 for te st

if sn==0; sv=0; else sv=1; end%set saving data to null if subject number =0
zbr=1; %activate zebar/optotune
% rc00=[0 0; 0 0]; %moving screen
% rc00=[40 0; 4 -7]; %
% rc00=[0 -8; 4 -7]; % updated april 20 2021, use this
rc00=[4 -7; 8 -6]; % updated June 210610, by fitting circle+centering it,looking thru system 210610

%  [right_Y right_X; left_Y Left_X]

sz = [1080, 1920]; % size of screen
screens = Screen('Screens');
screenNumber = max(screens);
load cal_val; %cf=[RB./RR LB./LR];
cf=ones(3,2);
%cf=[RB./RR LB./LR; RB./RG LB./LG; 1 1];
% cf=[0.6153 0.6153; 0.3982 0.3982; 1 1]; %gamma corrected purple
% background LATEST JAN 7
try;     
%     for k0=1:6; disp(opto(k0).control.focal_power); end
    for k0=1:6; disp(opto(k0).control.getFocalPower); end
    for k0=1:3; disp(zaber(k0).control.getposition); end
catch ERROR;
    for k0=1:3
        cls_opt; pause(3);
        
    end
    
    OPNT;
end
%%run up to here to START/OPEN the system


deg=-3;            
zaber(name_map('rotation')).move_deg(deg); %%-6400            
zaber(name_map('rotation')).control.getposition  


