%% SCRIPT FOR ANALYZING LCA DATA OF ALL PARTICIPANTS WHO PASSED SCREENING

clear;
dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';

%% MAKE FIGURE 7E

% subjNumAll = [1 3 5 9 10 16 17 18 20];
subjNumAll = [1 3 5 10 16 17 18 20];
% MEASURED DEFOCUS VALUES FOR EACH DISPLAY PRIMARY
defocusLCAmeasuredAll = [];
% FITTED LCA FUNCTION PARAMETERS
q1bestAll = []; 
q2bestAll = [];
q3bestAll = [];
% NEED THIS VARIABLE FOR ALIGNING LCA CURVES. THIS IS THE DEFOCUS AT THE 
% PEAK FOR THE GREEN PRIMARY (533nm) SPECIFIED BY THE FITTED LCA CURVE. 
% I USE IT TO SHIFT THE LCA CURVES VERTICALLY SO THEY ARE ALL ALIGNED AT 0, 
% THEN I SHIFT THEM AGAIN SO THEY ARE ALIGNED AT -2D. THEY WAY I DID IT IS 
% QUITE CONVOLUTED AND I SHOULD CLEAN THIS UP EVENTUALLY. BUT I'M SURE 
% IT'S CORRECT. 
D0all = []; 
% STORE BOOTSTRAPPED DEFOCUS VALUES FOR EACH DISPLAY PRIMARY
defocusLCAmeasuredBootsAll = [];

for i = 1:length(subjNumAll)
    % ANALYZING LCA FOR EACH PARTICIPANT
    [defocusLCAmeasured, q1best, q2best, q3best,defocusLCAmeasuredBoots,D0] = ARCacuAnalysisLCA(subjNumAll(i),0,100,dataPath);
    % STORING THREE DEFOCUS VALUES FOR RED, GREEN, AND BLUE
    defocusLCAmeasuredAll(i,:) = defocusLCAmeasured;
    % STORING PARAMETERS
    q1bestAll(i) = q1best;
    q2bestAll(i) = q2best;
    q3bestAll(i) = q3best;
    % STORING BOOTSTRAPS
    defocusLCAmeasuredBootsAll(:,:,i) = defocusLCAmeasuredBoots;
    D0all(i) = D0;
    display(['Finished subject ' num2str(i)]);
end

% NOT USING THE FOLLOWING 4 VARIABLES ANY LONGER (STANDARD LCA FUNCTION)
wave2fit = 460:5:875; % SUPPORT FOR LCA FUNCTION
standardLCAfit = humanWaveDefocus(wave2fit); % STANDARD LCA
val532 = humanWaveDefocus(532);
standardLCAfit = standardLCAfit-val532;

displacementLambda = -4:4; % JITTER DATA POINTS FOR EACH PARTICIPANT FOR VISIBILITY
wavePlot = 380:5:875; % SUPPORT OVER LCA FUNCTION
defocusAt875centered550 = []; % POORLY-NAMED VARIABLE

figure;
hold on;
% plot(wave2fit,standardLCAfit,'k-','LineWidth',1);
% PLOT CONTINUOUS LCA CURVES
for i = 1:size(defocusLCAmeasuredAll,1)
    set(gca,'ColorOrderIndex',i);
    % THE NAME OF THIS VARIABLE NO LONGER REFLECTS WHAT IT IS--I'LL COME
    % BACK TO RENAME THIS
    defocusAt875centered550(i) = humanWaveDefocusARC(533,wavePlot(end),subjNumAll(i)) + D0all(i)+defocusLCAmeasuredAll(i,2);
    plot(wavePlot,(-2-defocusAt875centered550(i))+humanWaveDefocusARC(533,wavePlot,subjNumAll(i))+D0all(i)+defocusLCAmeasuredAll(i,2),'-','LineWidth',1.5);
    defocus875to620(i,:) = humanWaveDefocusARC(533,wavePlot([48 100]),subjNumAll(i))+D0all(i)+defocusLCAmeasuredAll(i,2);
end
% PLOT ERROR BARS FOR THREE DATA POINTS (R, G, AND B) FOR EACH PARTICIPANT
for i = 1:size(defocusLCAmeasuredAll,1)
    defocusLCAmeasuredBootsTmp = squeeze(defocusLCAmeasuredBootsAll(:,:,i))-1*defocusLCAmeasuredAll(i,2);
    CIlca = quantile(defocusLCAmeasuredBootsTmp',[0.025 0.975]);
    defocusLCAmeasuredAll(i,:) = defocusLCAmeasuredAll(i,:);
    set(gca,'ColorOrderIndex',i);
    errorbar([616 533 468]-displacementLambda(i), ...
             -(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2))+(-2-defocusAt875centered550(i)), ...
             -(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2)-CIlca(2,:)), ...
             -(CIlca(1,:)-(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2))), ...
             '.','MarkerSize',10, ...
        'LineWidth',1.5,'MarkerFaceColor','w');    
end
% PLOT THREE DATA POINTS (R, G, AND B) FOR EACH PARTICIPANT
for i = 1:size(defocusLCAmeasuredAll,1)  
    set(gca,'ColorOrderIndex',i);
    plot([616 533 468]-displacementLambda(i),-(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2))+(-2-defocusAt875centered550(i)),'o','MarkerSize',10, ...
        'LineWidth',1.5,'MarkerFaceColor','w');
end
% axis square;
xlim([440 875]);
% ylim([-1.5 2]);
formatFigure('Wavelength (\lambda)','Defocus (D)');
legend('','','','','','','','','','','','','','','','','S1','S2','S3','S4','S5','S6','S7','S8','Location','SouthEast');

% figure;
% hold on;
% plot(ones([1 8]),defocus875to620(:,2)-defocus875to620(:,1),'ko','LineWidth',1);
% plot(1,humanWaveDefocus(875)-humanWaveDefocus(620),'k*','LineWidth',1);
% plot(1,mean(defocus875to620(:,2)-defocus875to620(:,1)),'kp','LineWidth',1);
% ylim([0 0.8]);
% ylabel('D_{875}-D_{620}');
% set(gca,'FontSize',15);
% set(gca,'XTick',[]);
