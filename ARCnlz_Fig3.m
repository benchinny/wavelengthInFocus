%% LOADED DATA GENERATED WITH FIRST PART OF ARCnlz_mainExpAvg

dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end

load([dataPath 'data' slash 'PresavedFigureData' slash 'allExp1DataRGB.mat']);

%% GENERATE FIGURE 3

% UNIQUE COLOR CONDITIONS IN TERMS OF PROPORTION OF MAX LUMINANCE
rgbLumNormCndUnq = [0.2500         0    1.0000; ...
                    0.5000         0    1.0000; ...
                    1.0000         0    1.0000; ...
                    1.0000         0    0.5000; ...
                    1.0000         0    0.2500; ...
                    0.2500    0.5000    1.0000; ...
                    0.5000    0.5000    1.0000; ...
                    1.0000    0.5000    1.0000; ...
                    1.0000    0.5000    0.5000; ...
                    1.0000    0.5000    0.2500; ...
                    1.0000    1.0000    1.0000];

nBoots = 1000;

defocusAt875objDistCenteredCell = {};
wvInFocusSortedCell = {};
optDistUnq = [1.5 2.5 3.5]; % UNIQUE OPTICAL DISTANCES
subjSymbols = 'sd+x^><v';

for k = 1:length(optDistUnq)
    for i = 1:size(rgbLumNormCndUnq,1)
        % GET ALL CELL INDICES WITH UNIQUE COLOR AND DISTANCE
        indRgbLumNormCnd = find(abs(rgbLumNormCndAll(:,1)-rgbLumNormCndUnq(i,1))<0.001 & ...
                                abs(rgbLumNormCndAll(:,2)-rgbLumNormCndUnq(i,2))<0.001 & ...
                                abs(rgbLumNormCndAll(:,3)-rgbLumNormCndUnq(i,3))<0.001 & ...
                                abs(optDistCndAll+optDistUnq(k))<0.001);
        defocusAt875objDistCentered = [];
        wvInFocusSorted = [];
        subjNumTagInd = subjNumTag(indRgbLumNormCnd);
        subjNumTagTmp = [];
        for j = 1:length(indRgbLumNormCnd)
            % DEFOCUS AT 875NM RELATIVE TO OBJECT DISTANCE
            defocusAt875objDistCentered = [defocusAt875objDistCentered; defocusAt875cellAll{indRgbLumNormCnd(j)}-optDistCndAll(indRgbLumNormCnd(j))];
            % SORTING WAVELENGTH-IN-FOCUS VALUES FROM PRE-SAVED DATA
            wvInFocusSorted = [wvInFocusSorted; wvInFocusCellAll{indRgbLumNormCnd(j)}];
            % SORTING DATA BY PARTICIPANT
            subjNumTagTmp = [subjNumTagTmp; subjNumTagInd(j).*ones(size(defocusAt875cellAll{indRgbLumNormCnd(j)}))];
        end
        subjNumTagUnq = unique(subjNumTagTmp);
        defocusAt875objDistCenteredAvg = [];
        wvInFocusAvg = [];
        % AVERAGING DEFOCUS AND WAVELENGTH-IN-FOCUS ACROSS TRIALS FOR EACH
        % PARTICIPANT
        for j = 1:length(subjNumTagUnq)
            defocusAt875objDistCenteredAvg(j) = mean(defocusAt875objDistCentered(subjNumTagTmp==subjNumTagUnq(j)));
            wvInFocusAvg(j) = mean(wvInFocusSorted(subjNumTagTmp==subjNumTagUnq(j)));
        end
        defocusAt875objDistCenteredCell{i} = defocusAt875objDistCenteredAvg;
        wvInFocusSortedCell{i,k} = wvInFocusAvg;
        bootInds = randsample(length(wvInFocusAvg),length(wvInFocusAvg)*nBoots,'true');
        bootIndsMatrix = reshape(bootInds,[length(wvInFocusAvg) nBoots]);
        wvInFocusBoots = mean(wvInFocusAvg(bootIndsMatrix),1);
        defocusAt875objDistCenteredBoots = mean(defocusAt875objDistCenteredAvg(bootIndsMatrix),1);
        wvInFocusCI(:,i) = quantile(wvInFocusBoots,[0.025 0.975]);
        defocusAt875CI(:,i) = quantile(defocusAt875objDistCenteredBoots,[0.025 0.975]);
        standardError95defocus(i) = 1.96*std(defocusAt875objDistCenteredAvg)/sqrt(length(defocusAt875objDistCenteredAvg));
        standardError95wvInFocus(i,k) = 1.96*std(wvInFocusAvg)/sqrt(length(wvInFocusAvg));
    end
end

distSymbols = 'sod';

figure;
set(gcf,'Position',[353 426 988 420]);
subplot(1,2,1);
hold on;
for i = 1:5 % size(rgbLumNormCndUnq,1)
    for k = 1:3
        plot(i,mean(wvInFocusSortedCell{i,k}),['k' distSymbols(k)],'MarkerFaceColor', ...
            rgbLumNormCndUnq(i,:),'MarkerSize',12,'LineWidth',1);
        wvInFocusSortedCellContent = wvInFocusSortedCell{i,k};
        % for j = 1:length(wvInFocusSortedCellContent)
        %     plot(i,wvInFocusSortedCellContent(j), ...
        %         subjSymbols(j),'MarkerSize',12,'Color',rgbLumNormCndUnq(i,:));    
        % end
    end
end
for i = 1:5 %size(rgbLumNormCndUnq,1)
    for k = 1:3
        CIrgb = quantile(wvInFocusSortedCell{i,k},[0.16 0.84]);
        % errorbar(i,mean(wvInFocusSortedCell{i}), ...
        %          mean(wvInFocusSortedCell{i})-wvInFocusCI(1,i), ...
        %          wvInFocusCI(2,i)-mean(wvInFocusSortedCell{i}),'ko','MarkerFaceColor', ...
        %         rgbLumNormCndUnq(i,:),'Color',rgbLumNormCndUnq(i,:),'MarkerSize',12,'LineWidth',1);
        errorbar(i,mean(wvInFocusSortedCell{i,k}), ...
                 standardError95wvInFocus(i,k),['k' distSymbols(k)],'MarkerFaceColor', ...
                rgbLumNormCndUnq(i,:),'Color',rgbLumNormCndUnq(i,:),'MarkerSize',12,'LineWidth',1);  
    end
end
axis square;
formatFigure('Red-blue ratio','Wavelength in focus (nm)','No green');
xlim([0.5 5.5]);
ylim([450 650]);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
subplot(1,2,2);
hold on;
for i = 1:5 % size(rgbLumNormCndUnq,1)
    for k = 1:3
        plot(i,mean(wvInFocusSortedCell{i+5,k}),['k' distSymbols(k)],'MarkerFaceColor', ...
            rgbLumNormCndUnq(i+5,:),'MarkerSize',12,'LineWidth',1);
        wvInFocusSortedCellContent = wvInFocusSortedCell{i+5,k};
        % for j = 1:length(subjSymbols)
        %     plot(i,wvInFocusSortedCellContent(j), ...
        %         subjSymbols(j),'MarkerSize',12,'Color',rgbLumNormCndUnq(i,:));        
        % end
    end
end
for i = 1:5 %size(rgbLumNormCndUnq,1)
    for k = 1:3
        CIrgb = quantile(wvInFocusSortedCell{i+5,k},[0.16 0.84]);
        % errorbar(i,mean(wvInFocusSortedCell{i+5}), ...
        %          mean(wvInFocusSortedCell{i+5})-wvInFocusCI(1,i+5), ...
        %          wvInFocusCI(2,i+5)-meclan(wvInFocusSortedCell{i+5}),'ko','MarkerFaceColor', ...
        %         rgbLumNormCndUnq(i+5,:),'Color',rgbLumNormCndUnq(i+5,:),'MarkerSize',12,'LineWidth',1);
        errorbar(i,mean(wvInFocusSortedCell{i+5,k}), ...
                 standardError95wvInFocus(i+5,k),['k' distSymbols(k)],'MarkerFaceColor', ...
                rgbLumNormCndUnq(i+5,:),'Color',rgbLumNormCndUnq(i+5,:),'MarkerSize',12,'LineWidth',1); 
    end
end
% plot([0.5 5.5],mean(defocusAt875objDistCenteredCell{end}).*[1 1],'k--','LineWidth',1);
axis square;
formatFigure('Red-blue ratio','Wavelength in focus (nm)','Some green');
xlim([0.5 5.5]);
ylim([450 650]);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
