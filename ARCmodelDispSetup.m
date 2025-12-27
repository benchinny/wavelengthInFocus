function d = ARCmodelDispSetup(calPath)

% sets up display struct in ISETBIO for ARchroma project
%
% only takes in path to display calibration files as input

% GAMMA EXPONENTS FOR EACH RGB CHANNEL OF DISPLAY
rGamma = 2.5;
gGamma = 2.7;
bGamma = 2.3;

% Setting up display properties
d = displayCreate('OLED-Samsung');
d = displaySet(d, 'name', 'my display');
% SIMULATED SCREEN DISTANCE--IT'S ARBITRARY SINCE THE BVAMS HAS
% NON-STANDARD OPTICAL PROPERTIES
d = displaySet(d,'ViewingDistance',1);
% 378 DOTS PER INCH YIELDS APPROXIMATELY 260 PIXELS PER VISUAL DEGREE,
% WHICH IS WHAT THE BVAMS HAS
d = displaySet(d,'dpi',378); 
% GET RID OF ALL UNNECESSARY FIELDS
d.dixel = [];
d.mainimage = [];

% LOAD BVAMS CALIBRATION DATA
load(fullfile(calPath,'redPrimaryJuly0624_initialPositionFocus3_100.mat'));
d.spd(:,1) = energy;
load(fullfile(calPath,'greenPrimaryJuly0624_initialPositionFocus3_100.mat'));
d.spd(:,2) = energy;
load(fullfile(calPath,'bluePrimaryJuly0624_initialPositionFocus3_100.mat'));
d.spd(:,3) = energy;

% APPLY OUR EMPIRICALLY DERIVED GAMMA FROM CALIBRATION MEASUREMENTS
d.gamma(:,1) = (linspace(0,1,1024)').^rGamma;
d.gamma(:,2) = (linspace(0,1,1024)').^gGamma;
d.gamma(:,3) = (linspace(0,1,1024)').^bGamma;

end