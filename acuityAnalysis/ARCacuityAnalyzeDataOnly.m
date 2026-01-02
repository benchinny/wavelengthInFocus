%%
function [unqFocDst,PC,PCci,dprime,dprimeCI,PCfit,dprimeFitAll,PCfitSupport,bestDist,bestDistCI,PCboots] = ARCacuityAnalyzeDataOnly(subjNum,bPLOT,dataPath)

% analyzes data from acuity experiment
%
% subjNum: subject ID. Valid numbers: 1, 3, 5, 10, 16, 17, 18, 20.
% bPLOT  : plot or not
% dataPath: path to data on local machine
%
% unqFocDst: unique stimulus optical distances
% PC: proportion correct
% PCci: confidence intervals on proportion correct
% dprime: d-prime values corresponding to PC
% dprimeCI: d-prime confidence intervals corresponding to PCci
% PCfit: spline fit to proportion correct
% dprimeFitAll: d-prime values corresponding to PCfit
% PCfitSupport: x-axis for PCfit
% bestDist: distance yielding best performance
% bestDistCI: confidence intervals on the best distance
% PCboots: bootstrapped proportion correct values

% WHERE THE ACUITY DATA RESIDE
dataDirectory = fullfile(dataPath,'data','psychophysicalData');
% SPECIFY 95% CONFIDENCE INTERVALS
confIntBounds = [0.025 0.975];

% FILENAMES FOR DATA
if subjNum==3
    filenames = {
                  % % S3 PURPLE (ACTUAL)
                  'S1013V19_AFC_RightACL0_2409201559.mat' ...
                  'S1013V19_AFC_RightACL0_2409201548.mat' ...
                  'S1013V19_AFC_RightACL0_2409201542.mat' ...
                  'S1013V19_AFC_RightACL0_2409201529.mat' ...   
                  };
elseif subjNum==1
    filenames = {
                  % % S1 PURPLE (ACTUAL)
                  'S1011V18_AFC_RightACL0_2411050943.mat' ...
                  'S1011V18_AFC_RightACL0_2411050928.mat' ...
                  };
elseif subjNum==10
    filenames = {
                  % S10 PURPLE
                  'S1020V10_AFC_RightACL0_2409111552.mat' ...
                  'S1020V10_AFC_RightACL0_2409111600.mat' ...
                  'S1020V10_AFC_RightACL0_2409111607.mat' ...
                  'S1020V10_AFC_RightACL0_2409111615.mat' ...
                  };
elseif subjNum==5
    filenames = {
                  % S10 PURPLE
                  'S1015V10_AFC_RightACL0_2410031638.mat' ...
                  'S1015V10_AFC_RightACL0_2410031631.mat' ...
                  'S1015V10_AFC_RightACL0_2410031615.mat' ...
                  'S1015V10_AFC_RightACL0_2410031608.mat' ...
                  };    
elseif subjNum==16
    filenames = {
                  % S16 PURPLE
                  'S1026V8_AFC_RightACL0_2409231008.mat' ...
                  'S1026V8_AFC_RightACL0_2409231000.mat' ...
                  'S1026V8_AFC_RightACL0_2409230953.mat' ...
                  'S1026V8_AFC_RightACL0_2409230946.mat' ...    
                  };    
elseif subjNum==17
    filenames = {
                  % S17 PURPLE
                  'S1027V9_AFC_RightACL0_2410081135.mat' ...
                  'S1027V9_AFC_RightACL0_2410081141.mat' ...
                  'S1027V9_AFC_RightACL0_2410081208.mat' ...
                  'S1027V9_AFC_RightACL0_2410081214.mat' ...    
                  }; 
elseif subjNum==18
    filenames = {
                  % S17 PURPLE
                  'S1028V9_AFC_RightACL0_2411051411.mat' ...
                  'S1028V9_AFC_RightACL0_2411051420.mat' ...
                  'S1028V9_AFC_RightACL0_2411151223.mat' ...
                  'S1028V9_AFC_RightACL0_2411151217.mat' ...    
                  }; 
elseif subjNum==20
    filenames = {
                  % S17 PURPLE
                  'S1030V8_AFC_RightACL0_2411151125.mat' ...
                  'S1030V8_AFC_RightACL0_2411151134.mat' ...
                  'S1030V8_AFC_RightACL0_2411151143.mat' ...
                  'S1030V8_AFC_RightACL0_2411151152.mat' ...    
                  };         
end

rgb = []; % COLOR CONDITIONS
meanFocstmOptDst = []; % DISTANCE OF ACCOMMODATIVE STIMULUS
focStmOptDstIncr = []; % DISTANCE OF ACUITY STIMULUS RELATIVE TO ACCOMMODATIVE STIMULUS
rspAcu = []; % RESPONSES: 1 FOR -15 DEG, 2 FOR 15 DEG
stimOrientation = []; % STIMULUS ORIENTATION: 1 FOR -15 DEG, 2 FOR 15 DEG

% SORTING DATA
for i = 1:length(filenames)
    % LOAD DATA
    load(fullfile(dataDirectory,filenames{i}));
    % STACK RELEVANT VARIABLES FROM EACH FILE
    rgb = [rgb; AFCp.rgb];
    meanFocstmOptDst = [meanFocstmOptDst; AFCp.meanFocstmOptDst];
    focStmOptDstIncr = [focStmOptDstIncr; AFCp.focStmOptDstIncr];
    rspAcu = [rspAcu; AFCp.rspAcu'];
    stimOrientation = [stimOrientation; AFCp.stimOrientation];
end

% UNIQUE COLOR CONDITIONS
rgbUnq = unique(rgb,'rows');

unqFocDst = unique(focStmOptDstIncr); % UNIQUE DISTANCES
% SCALE FACTOR FOR GOING FROM OPTOTUNE POWER CHANGE TO DEFOCUS AT EYE
scaleFac = 0.816;

for i = 1:length(unqFocDst) % LOOP OVER UNIQUE DISTANCE INCREMENTS OF ACUITY STIMULUS
     ind = focStmOptDstIncr==unqFocDst(i);
     % COMPUTE PROPORTION CORRECT
     PC(i) = sum(rspAcu(ind)==stimOrientation(ind))./sum(ind);
     % COMPUTE 95% CONFIDENCE INTERVALS
     PCci(:,i) = binoinv(confIntBounds,sum(ind),PC(i))./sum(ind);
end

% FITTING SPLINE
PCfitSupport = min(unqFocDst.*scaleFac):0.01:max(unqFocDst.*scaleFac);
PCfit = spline(unqFocDst.*scaleFac,PC,PCfitSupport);
% EXTRACT BEST DISTANCE
[~,indBest] = max(PCfit);
bestDist = PCfitSupport(indBest);

% BOOTSTRAPPING
nBoots = 1000; % HOW MANY BOOTSTRAPS
if nBoots>0
    for j = 1:nBoots % BOOTSTRAPPING
        for i = 1:length(unqFocDst)
            indAnalysis = focStmOptDstIncr==unqFocDst(i);
            % RESAMPLE INDICES (USING FIND TO CONVERT FROM BOOLEAN TO
            % NUMBER)
            indBoots = randsample(find(indAnalysis),sum(indAnalysis),'true');
            % COMPUTING PROPORTION CORRECT
            PCboots(i,j) = sum(rspAcu(indBoots)==stimOrientation(indBoots))./length(indBoots);
        end
        % FITTING SPLINE
        PCfitSupportBoots = min(unqFocDst.*scaleFac):0.01:max(unqFocDst.*scaleFac);
        PCfitBoots = spline(unqFocDst.*scaleFac,PCboots(:,j),PCfitSupportBoots);
        % FIND DISTANCE OF BEST ACUITY FOR EACH BOOTSTRAP
        [~,indBoots] = max(PCfitBoots);
        defocusMeasuredBoots(j) = PCfitSupportBoots(indBoots);
    end
end
% 95% CONFIDENCE INTERVALS ON BEST DISTANCE
bestDistCI = quantile(defocusMeasuredBoots,confIntBounds);

% FIRST, CONVERT PROPORTION CORRECT FITS TO D-PRIME
epsilonPC = 0.99; % ARTIFICIAL CAP TO PREVENT INFINITE D-PRIMES
PCfitDP = PCfit; % RENAME TO MAKE IT CLEAR WE'RE FITTING D-PRIME
% CAP PERFORMANCE TO PREVENT INFINITE D-PRIMES
PCfitDP(PCfitDP>epsilonPC) = epsilonPC;
% EQUATION FOR D-PRIME
dprimeFitAll = 2*norminv(PCfitDP);

% NOW, CALCULATE D-PRIME FOR ACTUAL DATA POINTS
PCforDP = PC;
% CAP PERFORMANCE TO PREVENT INFINITE D-PRIMES
PCforDP(PCforDP>epsilonPC) = epsilonPC;
% APPLY FORMULA FOR D-PRIME
dprime = 2*norminv(PCforDP);

% NOW, CALCULATE D-PRIME FOR CONFIDENCE INTERVALS
PCciForDP = PCci;
% CAP PERFORMANCE TO PREVENT INFINITE D-PRIMES
PCciForDP(PCciForDP>epsilonPC) = epsilonPC;
% APPLY FORMULA FOR D-PRIME
dprimeCI = 2*norminv(PCciForDP);

if bPLOT
    figure;
    hold on;
    % plot(unqFocDst.*scaleFac,PC,'o-','Color',rgbUnq,'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
    plot(PCfitSupport,PCfit,'-','Color',rgbUnq,'LineWidth',1.5);
    errorbar(unqFocDst.*scaleFac,PC,PC-PCci(1,:),PCci(2,:)-PC,'o','Color',rgbUnq,'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
    axis square;
    ylim([0.4 1]);
    formatFigure('Relative optical distance (D)','Proportion Correct',['Subject ' num2str(subjNum)]);
    
    figure;
    hold on;
    % plot(unqFocDst.*scaleFac,PC,'o-','Color',rgbUnq,'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
    plot(PCfitSupport,dprimeFitAll,'-','Color',rgbUnq,'LineWidth',1.5);
    errorbar(unqFocDst.*scaleFac,dprime,dprime-dprimeCI(1,:),dprimeCI(2,:)-dprime,'o','Color',rgbUnq,'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
    axis square;
    formatFigure('Relative optical distance (D)','d''',['Subject ' num2str(subjNum)]);
end

end
