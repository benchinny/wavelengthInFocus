function [q1bestAll, q2bestAll, q3bestAll] = ARCnlzLCAall(dataPath,bPLOT)

% FOR ANALYZING LCA DATA OF ALL PARTICIPANTS WHO PASSED SCREENING

% dataPath = 'C:\Users\bmccis\OneDrive - rit.edu\Documents\wavelengthInFocusData\';
%
% q1best: best fitting LCA curve parameter q1 (c in paper)
% q2best: best fitting LCA curve parameter q2 (a in paper)
% q3best: best fitting LCA curve parameter q3 (b in paper)

%% MAKE FIGURE 7E

% LIST OF SUBJECTS TO ANALYZE
subjNumAll = [1 3 5 10 16 17 18 20];
% MEASURED DEFOCUS VALUES FOR EACH DISPLAY PRIMARY
defocusLCAmeasuredAll = [];
% FITTED LCA FUNCTION PARAMETERS
q1bestAll = []; 
q2bestAll = [];
q3bestAll = [];

% STORE BOOTSTRAPPED DEFOCUS VALUES FOR EACH DISPLAY PRIMARY
defocusLCAmeasuredBootsAll = [];

for i = 1:length(subjNumAll)
    % ANALYZING LCA FOR EACH PARTICIPANT
    [defocusLCAmeasured, q1best, q2best, q3best,defocusLCAmeasuredBoots,~] = ARCnlzLCA(subjNumAll(i),0,500,dataPath);
    % STORING THREE DEFOCUS VALUES FOR RED, GREEN, AND BLUE
    defocusLCAmeasuredAll(i,:) = defocusLCAmeasured;
    % STORING PARAMETERS
    q1bestAll(i) = q1best;
    q2bestAll(i) = q2best;
    q3bestAll(i) = q3best;
    % STORING BOOTSTRAPS
    defocusLCAmeasuredBootsAll(:,:,i) = defocusLCAmeasuredBoots;
    display(['Finished subject ' num2str(i)]);
end

displacementLambda = -4:4; % JITTER DATA POINTS FOR EACH PARTICIPANT FOR VISIBILITY
wavePlot = 380:5:875; % SUPPORT OVER LCA FUNCTION

if bPLOT
    figure;
    hold on;
    % PLOT CONTINUOUS LCA CURVES
    amount2shift = [];
    for i = 1:size(defocusLCAmeasuredAll,1)
        set(gca,'ColorOrderIndex',i);
        % CALCULATE ORIGINAL DEFOCUS AT 875NM FOR FUNCTION FITS TO SUBJECT
        % DATA
        defocusAt875orig = humanWaveDefocusParameterized(875,q1bestAll(i),q2bestAll(i),q3bestAll(i));
        % CALCULATE SHIFT REQUIRED TO ANCHOR 875NM TO -2D
        amount2shift(i) = -2-defocusAt875orig;
        % PLOT WITH SHIFT INCLUDED
        plot(wavePlot,amount2shift(i)+humanWaveDefocusParameterized(wavePlot,q1bestAll(i),q2bestAll(i),q3bestAll(i)),'-','LineWidth',1.5);
    end
    % PLOT ERROR BARS FOR THREE DATA POINTS (R, G, AND B) FOR EACH PARTICIPANT
    for i = 1:size(defocusLCAmeasuredAll,1)
        % CALCULATE CONFIDENCE INTERVALS 
        CIlca = quantile(squeeze(defocusLCAmeasuredBootsAll(:,:,i))',[0.025 0.975]);
        set(gca,'ColorOrderIndex',i);
        % PLOT ERROR BARS ACCOUNTING FOR THE FACT THAT THE LCA FUNCTION IS
        % SIGN-FLIPPED RELATIVE TO THE DISTANCES OF BEST FOCUS FOR EACH
        % COLOR. 
        errorbar([616 533 468]-displacementLambda(i), ...
                 -defocusLCAmeasuredAll(i,:)+amount2shift(i), ...
                 -(defocusLCAmeasuredAll(i,:))+CIlca(2,:), ...
                 -(CIlca(1,:))+defocusLCAmeasuredAll(i,:), ...
                 '.','MarkerSize',10, ...
            'LineWidth',1.5,'MarkerFaceColor','w');    
    end
    % PLOT THREE DATA POINTS (R, G, AND B) FOR EACH PARTICIPANT
    for i = 1:size(defocusLCAmeasuredAll,1)  
        set(gca,'ColorOrderIndex',i);
        plot([616 533 468]-displacementLambda(i),-(defocusLCAmeasuredAll(i,:))+amount2shift(i),'o','MarkerSize',10, ...
            'LineWidth',1.5,'MarkerFaceColor','w');
    end
    % axis square;
    xlim([440 875]);
    ylim([-5 -1]);
    formatFigure('Wavelength (\lambda)','Defocus (D)');
    legend('','','','','','','','','','','','','','','','','S1','S2','S3','S4','S5','S6','S7','S8','Location','SouthEast');
end

% UNCOMMENT LINE BELOW TO SAVE
save(fullfile(dataPath,'data','PresavedFigureData','LCAparams.mat'),'q1bestAll','q2bestAll','q3bestAll');
