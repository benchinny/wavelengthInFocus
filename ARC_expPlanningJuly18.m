%%

rgbConditions = [0.555 0.320 1.00; ...
                 0.416 0.320 1.00; ...
                 0.312 0.320 1.00; ...
                 0.555 0.320 0.73; ...
                 0.555 0.320 0.533; ...
                 0.555 0.000 1.00; ...
                 0.416 0.000 1.00; ...
                 0.312 0.000 1.00; ...
                 0.555 0.000 0.73; ...
                 0.555 0.000 0.533; ...
                 ];

load('/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Meetings/MeetingJuly18/wvInFocusComparisonPermute2.mat');

reorderGreen = [5 4 1 2 3];

reorderNoGreen = [10 9 6 7 8];

figure;
set(gcf,'Position',[338 156 1087 788]);
subplot(2,2,1);
plot(1:length(reorderGreen),wvInFocusXCall(reorderGreen,:),'ko','MarkerSize',12, ...
     'MarkerFaceColor','w');
axis square;
set(gca,'FontSize',15);
xlabel('Condition');
ylabel('Wavelength in focus');
title('XC with green');
xlim([0.5 5.5]);
ylim([470 max(wvInFocusXCall(:))+10]);
set(gca,'XTickLabel',{});

subplot(2,2,2);
plot(1:length(reorderNoGreen),wvInFocusXCall(reorderNoGreen,:),'ko','MarkerSize',12, ...
     'MarkerFaceColor','w');
axis square;
set(gca,'FontSize',15);
xlabel('Condition');
ylabel('Wavelength in focus');
title('XC without green');
xlim([0.5 5.5]);
ylim([470 max(wvInFocusXCall(:))+10]);set(gca,'XTickLabel',{});

subplot(2,2,3);
plot(1:length(reorderGreen),wvInFocusSTall(reorderGreen,:),'ko','MarkerSize',12, ...
     'MarkerFaceColor','w');
axis square;
set(gca,'FontSize',15);
xlabel('Condition');
ylabel('Wavelength in focus');
title('Strehl with green');
xlim([0.5 5.5]);
ylim([470 max(wvInFocusSTall(:))+10]);
set(gca,'XTickLabel',{});

subplot(2,2,4);
plot(1:length(reorderNoGreen),wvInFocusSTall(reorderNoGreen,:),'ko','MarkerSize',12, ...
     'MarkerFaceColor','w');
axis square;
set(gca,'FontSize',15);
xlabel('Condition');
ylabel('Wavelength in focus');
title('Strehl without green');
xlim([0.5 5.5]);
ylim([470 max(wvInFocusSTall(:))+10]);
set(gca,'XTickLabel',{});

figure;
set(gcf,'Position',[338 156 1087 788]);
subplot(2,2,1);
plot(1:length(reorderGreen),humanWaveDefocus(wvInFocusXCall(reorderGreen,:)),'ko','MarkerSize',12, ...
     'MarkerFaceColor','w');
axis square;
set(gca,'FontSize',15);
xlabel('Condition');
ylabel('Relative focus (D)');
title('XC with green');
xlim([0.5 5.5]);
ylim([-0.8 0.2]);
set(gca,'XTickLabel',{});

subplot(2,2,2);
plot(1:length(reorderNoGreen),humanWaveDefocus(wvInFocusXCall(reorderNoGreen,:)),'ko','MarkerSize',12, ...
     'MarkerFaceColor','w');
axis square;
set(gca,'FontSize',15);
xlabel('Condition');
ylabel('Relative focus (D)');
title('XC without green');
xlim([0.5 5.5]);
ylim([-0.8 0.2]);
set(gca,'XTickLabel',{});

subplot(2,2,3);
plot(1:length(reorderGreen),humanWaveDefocus(wvInFocusSTall(reorderGreen,:)),'ko','MarkerSize',12, ...
     'MarkerFaceColor','w');
axis square;
set(gca,'FontSize',15);
xlabel('Condition');
ylabel('Relative focus (D)');
title('Strehl with green');
xlim([0.5 5.5]);
ylim([min(humanWaveDefocus(wvInFocusSTall(:)))-0.05 max(humanWaveDefocus(wvInFocusSTall(:)))+0.05]);
set(gca,'XTickLabel',{});

subplot(2,2,4);
plot(1:length(reorderNoGreen),humanWaveDefocus(wvInFocusSTall(reorderNoGreen,:)),'ko','MarkerSize',12, ...
     'MarkerFaceColor','w');
axis square;
set(gca,'FontSize',15);
xlabel('Condition');
ylabel('Relative focus (D)');
title('Strehl without green');
xlim([0.5 5.5]);
ylim([min(humanWaveDefocus(wvInFocusSTall(:)))-0.05 max(humanWaveDefocus(wvInFocusSTall(:)))+0.05]);
set(gca,'XTickLabel',{});

