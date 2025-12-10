%% GENERATE FIGURE 3

% LOADED DATA GENERATED WITH FIRST PART OF ARCnlz_mainExpAvg
dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
if ispc
    slash = '\';
else
    slash = '/';
end

load([dataPath 'data' slash 'PresavedFigureData' slash 'allExp1DataRGB.mat']);

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

wvInFocusSortedCell = {}; % MAIN THING TO PLOT
optDistUnq = [1.5 2.5 3.5]; % UNIQUE OPTICAL DISTANCES
subjSymbols = 'sd+x^><v';

for k = 1:length(optDistUnq) % LOOP OVER OPTICAL DISTANCES
    for i = 1:size(rgbLumNormCndUnq,1) % LOOP OVER COLOR CONDITION
        % GET ALL CELL INDICES WITH UNIQUE COLOR AND DISTANCE
        indRgbLumNormCnd = find(abs(rgbLumNormCndAll(:,1)-rgbLumNormCndUnq(i,1))<0.001 & ...
                                abs(rgbLumNormCndAll(:,2)-rgbLumNormCndUnq(i,2))<0.001 & ...
                                abs(rgbLumNormCndAll(:,3)-rgbLumNormCndUnq(i,3))<0.001 & ...
                                abs(optDistCndAll+optDistUnq(k))<0.001);
        % THESE VALUES ARE ALL DESCRIBED AFTER THE NEXT 'FOR' LOOP
        wvInFocusSorted = [];
        subjNumTagInd = subjNumTag(indRgbLumNormCnd);
        subjNumTagTmp = [];
        for j = 1:length(indRgbLumNormCnd)
            % SORTING WAVELENGTH-IN-FOCUS VALUES FROM PRE-SAVED DATA
            wvInFocusSorted = [wvInFocusSorted; wvInFocusCellAll{indRgbLumNormCnd(j)}];
            % SORTING DATA BY PARTICIPANT
            subjNumTagTmp = [subjNumTagTmp; subjNumTagInd(j).*ones(size(defocusAt875cellAll{indRgbLumNormCnd(j)}))];
        end
        subjNumTagUnq = unique(subjNumTagTmp); % UNIQUE SUBJECT TAGS
        wvInFocusAvg = [];
        % WAVELENGTH-IN-FOCUS ACROSS TRIALS FOR EACH PARTICIPANT
        for j = 1:length(subjNumTagUnq) % LOOP OVER SUBJECTS
            % AVERAGE WAVELENGTH-IN-FOCUS
            wvInFocusAvg(j) = mean(wvInFocusSorted(subjNumTagTmp==subjNumTagUnq(j)));
        end
        % STORE IN A CELL
        wvInFocusSortedCell{i,k} = wvInFocusAvg;
        % STD ERROR FOR PLOTTING
        standardError95wvInFocus(i,k) = 1.96*std(wvInFocusAvg)/sqrt(length(wvInFocusAvg));
    end
end

distSymbols = 'sod';

figure;
set(gcf,'Position',[353 426 988 420]);
subplot(1,2,1);
hold on;
% PLOT DATA POINTS
for i = 1:5 % FOR 'NO GREEN' CONDITIONS
    for k = 1:3 % FOR ALL OPTICAL DISTANCES
        plot(i,mean(wvInFocusSortedCell{i,k}),['k' distSymbols(k)],'MarkerFaceColor', ...
            rgbLumNormCndUnq(i,:),'MarkerSize',12,'LineWidth',1);
        wvInFocusSortedCellContent = wvInFocusSortedCell{i,k};
    end
end
% PLOT ERROR BARS
for i = 1:5 % FOR 'NO GREEN' CONDITIONS
    for k = 1:3 % FOR ALL OPTICAL DISTANCES
        CIrgb = quantile(wvInFocusSortedCell{i,k},[0.16 0.84]);
        errorbar(i,mean(wvInFocusSortedCell{i,k}), ...
                 standardError95wvInFocus(i,k),['k' distSymbols(k)],'MarkerFaceColor', ...
                rgbLumNormCndUnq(i,:),'Color',rgbLumNormCndUnq(i,:),'MarkerSize',12,'LineWidth',1);  
    end
end
axis square;
formatFigure('Red-blue ratio','Wavelength in focus (nm)','No green');
xlim([0.5 5.5]);
ylim([450 680]);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});
subplot(1,2,2);
hold on;
% PLOT DATA POINTS
for i = 1:5 % FOR 'SOME GREEN' CONDITIONS
    for k = 1:3 % FOR ALL OPTICAL DISTANCES
        plot(i,mean(wvInFocusSortedCell{i+5,k}),['k' distSymbols(k)],'MarkerFaceColor', ...
            rgbLumNormCndUnq(i+5,:),'MarkerSize',12,'LineWidth',1);
        wvInFocusSortedCellContent = wvInFocusSortedCell{i+5,k};
    end
end
% PLOT ERROR BARS
for i = 1:5 % FOR 'SOME GREEN' CONDITIONS
    for k = 1:3 % FOR ALL OPTICAL DISTANCES
        CIrgb = quantile(wvInFocusSortedCell{i+5,k},[0.16 0.84]);
        errorbar(i,mean(wvInFocusSortedCell{i+5,k}), ...
                 standardError95wvInFocus(i+5,k),['k' distSymbols(k)],'MarkerFaceColor', ...
                rgbLumNormCndUnq(i+5,:),'Color',rgbLumNormCndUnq(i+5,:),'MarkerSize',12,'LineWidth',1); 
    end
end
axis square;
formatFigure('Red-blue ratio','Wavelength in focus (nm)','Some green');
xlim([0.5 5.5]);
ylim([450 680]);
set(gca,'XTickLabel',{'0.25' '0.50' '1.00' '2.00' '4.00'});