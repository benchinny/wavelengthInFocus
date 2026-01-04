%% Script for recreating plots from Figure 7B in manuscript

% WHERE DATA IS LOCATED
dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
% SET SUBJECT NUMBER--USE EITHER 3 OR 10 (S2 OR S4) TO GET PLOTS SHOWN IN
% FIGURE 7B IN MANUSCRIPT
subjNum = 3;
% RGB VALUES OF EXAMPLE STIMULI TO PLOT
redColor = [0.569 0.334 0.547];
blueColor = [0.327 0.334 1];
% TRANSFORM COLOR FOR PLOTTING
rgbLumNormRed = [];
rgbLumNormRed(:,1) = (redColor(:,1).^2.5)./0.2442;
rgbLumNormRed(:,2) = (redColor(:,2).^2.7)./0.1037;
rgbLumNormRed(:,3) = (redColor(:,3).^2.3)./1;
rgbLumNormRed(rgbLumNormRed>1) = 1;
% TRANSFORM COLOR FOR PLOTTING
rgbLumNormBlue = [];
rgbLumNormBlue(:,1) = (blueColor(:,1).^2.5)./0.2442;
rgbLumNormBlue(:,2) = (blueColor(:,2).^2.7)./0.1037;
rgbLumNormBlue(:,3) = (blueColor(:,3).^2.3)./1;
rgbLumNormBlue(rgbLumNormBlue>1) = 1;

% GRAB RAW DATA FOR REDDEST TRIALS
[defocus875rawCellRed, timeStampRawCellRed] = ARCmakeFig7Bhelper(subjNum,redColor,2.5,dataPath);

% GRAB RAW DATA FOR BLUEST TRIALS
[defocus875rawCellBlue, timeStampRawCellBlue] = ARCmakeFig7Bhelper(subjNum,blueColor,2.5,dataPath);

% FOR CATCHING TRIALS WHERE THE WAVEFRONT SENSOR RAN A BIT OVER TIME
maxLength = 90;
trialMeanRed = []; % FOR STORING MEANS
trialMeanBlue = []; % FOR STORING MEANS

figure;
hold on;
for i = 1:length(defocus875rawCellRed)
    timestampTmp = timeStampRawCellRed{i};
    defocus875tmp = -defocus875rawCellRed{i};
    if length(timestampTmp)>90
        timestampTmp = timestampTmp(1:90);
        defocus875tmp = defocus875tmp(1:90);
    end
    plot(timestampTmp,defocus875tmp,'Color',rgbLumNormRed,'LineWidth',1);
    trialMeanRed(i) = mean(defocus875tmp);
end
plot(xlim,mean(trialMeanRed).*[1 1],'Color',rgbLumNormRed,'LineWidth',1);
for i = 1:length(defocus875rawCellBlue)
    timestampTmp = timeStampRawCellBlue{i};
    defocus875tmp = -defocus875rawCellBlue{i};
    if length(timestampTmp)>90
        timestampTmp = timestampTmp(1:90);
        defocus875tmp = defocus875tmp(1:90);
    end    
    plot(timestampTmp,defocus875tmp,'Color',rgbLumNormBlue,'LineWidth',1);
    trialMeanBlue(i) = mean(defocus875tmp);
end
plot(xlim,mean(trialMeanBlue).*[1 1],'Color',rgbLumNormBlue,'LineWidth',1);
axis square;
ylim([-2.2 -0.7]);
xlim([0 3.25]);
set(gca,'YTick',[-2.0 -1.5 -1.0]);
formatFigure('Time (s)','Defocus (D)');