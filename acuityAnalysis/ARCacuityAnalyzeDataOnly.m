%%
function [unqFocDst,PC,PCci,dprime,dprimeCI,PCfit,dprimeFitAll,PCfitSupport,bestDist,bestDistCI,PCboots] = ARCacuAnalysisSubjective(subjNum,bPLOT,dataPath)

if ispc
    slash = '\';
else
    slash = '/';
end
dataDirectory = [dataPath 'data' slash 'ARC' slash];

if subjNum==3
    filenames = {
                  % % S3 PURPLE (ACTUAL)
                  [dataDirectory 'S1013V19_AFC_RightACL0_2409201559.mat'] ...
                  [dataDirectory 'S1013V19_AFC_RightACL0_2409201548.mat'] ...
                  [dataDirectory 'S1013V19_AFC_RightACL0_2409201542.mat'] ...
                  [dataDirectory 'S1013V19_AFC_RightACL0_2409201529.mat'] ...   
                  };
elseif subjNum==1
    filenames = {
                  % % S1 PURPLE (ACTUAL)
                  [dataDirectory 'S1011V18_AFC_RightACL0_2411050943.mat'] ...
                  [dataDirectory 'S1011V18_AFC_RightACL0_2411050928.mat'] ...
                  };
elseif subjNum==9
    filenames = {
                  % S10 PURPLE
                  [dataDirectory 'S1019V9_AFC_RightACL0_2410041308.mat'] ...
                  [dataDirectory 'S1019V9_AFC_RightACL0_2410041301.mat'] ...
                  [dataDirectory 'S1019V9_AFC_RightACL0_2410041254.mat'] ...
                  [dataDirectory 'S1019V9_AFC_RightACL0_2410041248.mat'] ...
                  };    
elseif subjNum==10
    filenames = {
                  % S10 PURPLE
                  [dataDirectory 'S1020V10_AFC_RightACL0_2409111552.mat'] ...
                  [dataDirectory 'S1020V10_AFC_RightACL0_2409111600.mat'] ...
                  [dataDirectory 'S1020V10_AFC_RightACL0_2409111607.mat'] ...
                  [dataDirectory 'S1020V10_AFC_RightACL0_2409111615.mat'] ...
                  };
elseif subjNum==5
    filenames = {
                  % S10 PURPLE
                  [dataDirectory 'S1015V10_AFC_RightACL0_2410031638.mat'] ...
                  [dataDirectory 'S1015V10_AFC_RightACL0_2410031631.mat'] ...
                  [dataDirectory 'S1015V10_AFC_RightACL0_2410031615.mat'] ...
                  [dataDirectory 'S1015V10_AFC_RightACL0_2410031608.mat'] ...
                  };    
elseif subjNum==16
    filenames = {
                  % S16 PURPLE
                  [dataDirectory 'S1026V8_AFC_RightACL0_2409231008.mat'] ...
                  [dataDirectory 'S1026V8_AFC_RightACL0_2409231000.mat'] ...
                  [dataDirectory 'S1026V8_AFC_RightACL0_2409230953.mat'] ...
                  [dataDirectory 'S1026V8_AFC_RightACL0_2409230946.mat'] ...    
                  };    
elseif subjNum==17
    filenames = {
                  % S17 PURPLE
                  [dataDirectory 'S1027V9_AFC_RightACL0_2410081135.mat'] ...
                  [dataDirectory 'S1027V9_AFC_RightACL0_2410081141.mat'] ...
                  [dataDirectory 'S1027V9_AFC_RightACL0_2410081208.mat'] ...
                  [dataDirectory 'S1027V9_AFC_RightACL0_2410081214.mat'] ...    
                  }; 
elseif subjNum==18
    filenames = {
                  % S17 PURPLE
                  [dataDirectory 'S1028V9_AFC_RightACL0_2411051411.mat'] ...
                  [dataDirectory 'S1028V9_AFC_RightACL0_2411051420.mat'] ...
                  [dataDirectory 'S1028V9_AFC_RightACL0_2411151223.mat'] ...
                  [dataDirectory 'S1028V9_AFC_RightACL0_2411151217.mat'] ...    
                  }; 
elseif subjNum==20
    filenames = {
                  % S17 PURPLE
                  [dataDirectory 'S1030V8_AFC_RightACL0_2411151125.mat'] ...
                  [dataDirectory 'S1030V8_AFC_RightACL0_2411151134.mat'] ...
                  [dataDirectory 'S1030V8_AFC_RightACL0_2411151143.mat'] ...
                  [dataDirectory 'S1030V8_AFC_RightACL0_2411151152.mat'] ...    
                  };         
end

rgb = [];
meanFocstmOptDst = [];
focStmOptDstIncr = [];
rspAcu = [];
stimOrientation = [];

% SORTING DATA
for i = 1:length(filenames)
    load(filenames{i});
    rgb = [rgb; AFCp.rgb];
    meanFocstmOptDst = [meanFocstmOptDst; AFCp.meanFocstmOptDst];
    focStmOptDstIncr = [focStmOptDstIncr; AFCp.focStmOptDstIncr];
    rspAcu = [rspAcu; AFCp.rspAcu'];
    stimOrientation = [stimOrientation; AFCp.stimOrientation];
end

if size(unique(rgb,'rows'),1)>1
   error('ARCacuAnalysisPilot: code does not handle multiple colors yet!');
else
   rgbUnq = unique(rgb,'rows');
end

unqFocDst = unique(focStmOptDstIncr);
scaleFac = 0.816;

for i = 1:length(unqFocDst)
%    PC(i) = sum(rspAcu(focStmOptDstIncr==unqFocDst(i) & meanFocstmOptDst==meanFocInt)==stimOrientation(focStmOptDstIncr==unqFocDst(i) & meanFocstmOptDst==meanFocInt))./sum(focStmOptDstIncr==unqFocDst(i) & meanFocstmOptDst==meanFocInt); 
     PC(i) = sum(rspAcu(focStmOptDstIncr==unqFocDst(i))==stimOrientation(focStmOptDstIncr==unqFocDst(i)))./sum(focStmOptDstIncr==unqFocDst(i));
     PCci(:,i) = binoinv([0.025 0.975],sum(focStmOptDstIncr==unqFocDst(i)),PC(i))./sum(focStmOptDstIncr==unqFocDst(i));
end

% FITTING SPLINE
PCfitSupport = min(unqFocDst.*scaleFac):0.01:max(unqFocDst.*scaleFac);
PCfit = spline(unqFocDst.*scaleFac,PC,PCfitSupport);
[~,indBest] = max(PCfit);
bestDist = PCfitSupport(indBest);

nBoots = 1000;
if nBoots>0
    for j = 1:nBoots % BOOTSTRAPPING
        for i = 1:length(unqFocDst)
            indAnalysis = focStmOptDstIncr==unqFocDst(i);
            indBoots = randsample(find(indAnalysis),sum(indAnalysis),'true');
            PCboots(i,j) = sum(rspAcu(indBoots)==stimOrientation(indBoots))./length(indBoots);
        end
        PCfitSupportBoots = min(unqFocDst.*scaleFac):0.01:max(unqFocDst.*scaleFac);
        PCfitBoots = spline(unqFocDst.*scaleFac,PCboots(:,j),PCfitSupportBoots);
        [~,indBoots] = max(PCfitBoots);
        defocusMeasuredBoots(j) = PCfitSupportBoots(indBoots);
    end
end
bestDistCI = quantile(defocusMeasuredBoots,[0.025 0.975]);

% CONVERTING TO D-PRIME
epsilonPC = 0.99;
PCfitDP = PCfit;
PCfitDP(PCfitDP>epsilonPC) = epsilonPC;
dprimeFitAll = 2*norminv(PCfitDP);
PCforDP = PC;
PCforDP(PCforDP>epsilonPC) = epsilonPC;
dprime = 2*norminv(PCforDP);
PCciForDP = PCci;
PCciForDP(PCciForDP>epsilonPC) = epsilonPC;
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
