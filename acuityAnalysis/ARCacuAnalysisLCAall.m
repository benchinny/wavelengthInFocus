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
defocusAt875centered533 = []; % POORLY-NAMED VARIABLE

figure;
hold on;
% plot(wave2fit,standardLCAfit,'k-','LineWidth',1);
% PLOT CONTINUOUS LCA CURVES
for i = 1:size(defocusLCAmeasuredAll,1)
    set(gca,'ColorOrderIndex',i);
    % FIRST, WE COMPUTE THE DIFFERENCE BETWEEN THE
    % PEAK OF THE GREEN PRIMARY (533NM) AND 875NM (POSITIVE NUMBER). THEN
    % WE ADD D0all (THE VALUE AT 533NM FOR THE FIT LCA FUNCTION) TO GET THE
    % ACTUAL VALUE AT 875NM FOR THE FIT LCA FUNCTION. FINALLY, WE SUBTRACT
    % THE DEFOCUS CORRESPONDING TO THE EMPIRICAL MEASUREMENT AT 533NM
    % (MEASUREMENT VALUES ARE SIGN-FLIPPED RELATIVE TO THE LCA FUNCTION, SO
    % IT LOOKS LIKE WE'RE ADDING, BUT WE'RE ACTUALLY SUBTRACTING). SO THIS
    % IS THE DEFOCUS AT 875NM CENTERED ON 533NM.
    defocusAt875centered533(i) = humanWaveDefocusARC(533,wavePlot(end),subjNumAll(i)) + D0all(i)+defocusLCAmeasuredAll(i,2);
    % EXPLANATION OF LINE BELOW: TO GENERATE LCA CURVE ANCHORED AT 533,
    % CALL humanWaveDefocusARC(533,wavePlot,subjNumAll(i)). TO ANCHOR TO
    % SUBJECT'S ACTUAL DEFOCUS AT 533 AND THEN ANCHOR ENTIRE CURVE SUCH
    % THAT SUBJECT'S GREEN (533NM) DATAPOINT IS AT 0, CALL
    % +D0all(i)+defocusLCAmeasuredAll(i,2). FOR THE RESULTING FUNCTIONS,
    % ITS DEFOCUS AT 875NM WILL BE GIVEN BY THE VARIABLE
    % defocusAt875centered533. TO ANCHOR TO 2, SUBTRACT THAT VARIABLE, THEN
    % SUBTRACT 2.
    plot(wavePlot,(-2-defocusAt875centered533(i))+humanWaveDefocusARC(533,wavePlot,subjNumAll(i))+D0all(i)+defocusLCAmeasuredAll(i,2),'-','LineWidth',1.5);
end
% PLOT ERROR BARS FOR THREE DATA POINTS (R, G, AND B) FOR EACH PARTICIPANT
for i = 1:size(defocusLCAmeasuredAll,1)
    % SUBTRACT THE GREEN (533NM) DATAPOINT FROM ALL BOOTSTRAPS. 
    defocusLCAmeasuredBootsTmp = squeeze(defocusLCAmeasuredBootsAll(:,:,i))-1*defocusLCAmeasuredAll(i,2);
    CIlca = quantile(defocusLCAmeasuredBootsTmp',[0.025 0.975]);
    defocusLCAmeasuredAll(i,:) = defocusLCAmeasuredAll(i,:);
    set(gca,'ColorOrderIndex',i);
    % FIRST, ANCHOR THREE DATA POINTS SO THAT THE GREEN (533NM) ONE IS AT
    % 0. THIS WILL ALIGN ALL DATAPOINTS WITH THE LCA FUNCTION ANCHORED AT 
    % 0. THEN SUBTRACT 2 AND defocusAt875centered533 SO THAT THE 875NM
    % POINT IS AT -2 DIOPTERS. FOR NEGATIVE ERROR BAR LENGTH, ANCHOR
    % EMPIRICAL DEFOCUS MEASUREMENTS TO 0, THEN SUBTRACT THE UPPER CI. WE
    % DO THE UPPER CI BECAUSE IT WILL BECOME THE LOWER CI ONCE THE SIGN IS
    % FLIPPED. REPEAT FOR POSITIVE ERROR BAR LENGTH, BUT WITH EVERYTHING
    % FLIPPED. 
    errorbar([616 533 468]-displacementLambda(i), ...
             -(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2))+(-2-defocusAt875centered533(i)), ...
             -(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2)-CIlca(2,:)), ...
             -(CIlca(1,:)-(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2))), ...
             '.','MarkerSize',10, ...
        'LineWidth',1.5,'MarkerFaceColor','w');    
end
% PLOT THREE DATA POINTS (R, G, AND B) FOR EACH PARTICIPANT
for i = 1:size(defocusLCAmeasuredAll,1)  
    set(gca,'ColorOrderIndex',i);
    plot([616 533 468]-displacementLambda(i),-(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2))+(-2-defocusAt875centered533(i)),'o','MarkerSize',10, ...
        'LineWidth',1.5,'MarkerFaceColor','w');
end
% axis square;
xlim([440 875]);
% ylim([-1.5 2]);
formatFigure('Wavelength (\lambda)','Defocus (D)');
legend('','','','','','','','','','','','','','','','','S1','S2','S3','S4','S5','S6','S7','S8','Location','SouthEast');

