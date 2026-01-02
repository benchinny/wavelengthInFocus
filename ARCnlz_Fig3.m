%% MAKE FIGURE 3 IN THE PAPER

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
foldername = fullfile(dataPath,'data','PresavedFigureData');

% PRESAVED DATA:
load(fullfile(foldername,'allExp1DataRGB.mat'));

symbDist = 'sod'; % SYMBOLS FOR PLOTTING
% ORDER CONDITIONS FOR PLOTTING
conditionsOrderedNorm = [0.25 0.00 1.00; ...
                         0.50 0.00 1.00; ...
                         1.00 0.00 1.00; ...
                         1.00 0.00 0.50; ...
                         1.00 0.00 0.25; ...
                         0.25 0.50 1.00; ...
                         0.50 0.50 1.00; ...
                         1.00 0.50 1.00; ...
                         1.00 0.50 0.50; ...
                         1.00 0.50 0.25; ...
                         ];

subjNum = [1 3 5 10 16 17 18 20];

figure;
set(gcf,'Position',[353 426 988 420]);
% PLOT 'NO GREEN' CONDITIONS
subplot(1,2,1);
hold on;
for i = 1:size(wvMeanAll,2)
    wvMeanTmp = squeeze(wvMeanAll(:,i,:));
    for j = 1:5
        errorbar(j,mean(wvMeanTmp(j,:),2),1.96.*std(wvMeanTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
    end
end
axis square;
set(gca,'FontSize',15);
set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
set(gca,'Box','on');
xlim([0.5 5.5]);
ylim([450 680]);
xlabel('Red-blue ratio');
ylabel('Wavelength in focus (nm)');
% PLOT 'SOME GREEN' CONDITIONS
subplot(1,2,2);
hold on;
for i = 1:size(wvMeanAll,2)
    wvMeanTmp = squeeze(wvMeanAll(:,i,:));
    for j = 6:10
        errorbar(j-5,mean(wvMeanTmp(j,:),2),1.96.*std(wvMeanTmp(j,:)')./sqrt(8),['k' symbDist(i)],'MarkerSize',10,'MarkerFaceColor',conditionsOrderedNorm(j,:),'LineWidth',1);
    end
end
axis square;
xlim([0.5 5.5]);
ylim([450 680]);
set(gca,'Box','on');
set(gca,'FontSize',15);
set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
xlabel('Red-blue ratio');