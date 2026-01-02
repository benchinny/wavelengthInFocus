function [defocusLCAmeasured, q1best, q2best, q3best,defocusLCAmeasuredBoots,Dgreen,dprimeFitAll,PCfitSupport] = ARCnlzLCA(subj,bPLOT,nBoots,dataPath)

% analyzes data from the LCA experiment. Bootstraps data as well. 
%
% subj: subject numbers. Valid subj numbers: 1, 3, 5, 10, 16, 17, 18, 20.
% bPLOT: plot or not
% nBoots: number of bootstraps
% dataPath: data path
%
% defocusLCAmeasured: empirical defocus values for R, G, and B Gabors
% q1best: best fitting LCA curve parameter q1 (c in paper)
% q2best: best fitting LCA curve parameter q2 (a in paper)
% q3best: best fitting LCA curve parameter q3 (b in paper)
% defocusLCAmeasuredBoots: bootstrapped defocus values
% Dgreen: defocus value at green
% dprimeFitAll: d-prime values associated with fits
% PCfitSupportAll: support over which to fit proportion correct

rng(1); % INITIALIZE SAME RANDOM SEED

% DATA DIRECTORY
dataDirectory = fullfile(dataPath,'data','psychophysicalData');

% DATA TO LOAD
if subj==1
    filenames = {
                 'S1011V18_AFC_RightACL0_2410290922.mat' ...
                 'S1011V18_AFC_RightACL0_2410290936.mat' ...
                 'S1011V18_AFC_RightACL0_2410290949.mat' ...  
                 'S1011V18_AFC_RightACL0_2410291002.mat' ...    
                 'S1011V18_AFC_RightACL0_2410291014.mat' ...     
                 'S1011V18_AFC_RightACL0_2410291029.mat' ...
                 'S1011V18_AFC_RightACL0_2411050914.mat' ...
                 }; 
elseif subj==3
    filenames = {
                 'S1013V19_AFC_RightACL0_2409191523.mat' ...
                 'S1013V19_AFC_RightACL0_2409191540.mat' ...
                 'S1013V19_AFC_RightACL0_2409191547.mat' ...
                 'S1013V19_AFC_RightACL0_2409191554.mat' ...
                 'S1013V19_AFC_RightACL0_2409191609.mat' ... 
                 'S1013V19_AFC_RightACL0_2409191616.mat' ... 
                 'S1013V19_AFC_RightACL0_2409191631.mat' ...
                 'S1013V19_AFC_RightACL0_2409191640.mat' ... 
                 'S1013V19_AFC_RightACL0_2409191650.mat' ...     
                 'S1013V19_AFC_RightACL0_2409201521.mat' ...
                 'S1013V19_AFC_RightACL0_2409201515.mat' ... 
                 'S1013V19_AFC_RightACL0_2409201508.mat' ...                      
                 };            
elseif subj==10
    filenames = {
                 'S1020V11_AFC_RightACL0_2409171805.mat' ...
                 'S1020V11_AFC_RightACL0_2409171756.mat' ...
                 'S1020V11_AFC_RightACL0_2409171749.mat' ...
                 'S1020V11_AFC_RightACL0_2409171741.mat' ...
                 'S1020V11_AFC_RightACL0_2409171734.mat' ... 
                 'S1020V11_AFC_RightACL0_2409171728.mat' ... 
                 'S1020V11_AFC_RightACL0_2409171722.mat' ...
                 'S1020V11_AFC_RightACL0_2409171715.mat' ... 
                 'S1020V11_AFC_RightACL0_2409171707.mat' ... 
                 'S1020V11_AFC_RightACL0_2409171701.mat' ... 
                 'S1020V11_AFC_RightACL0_2409171654.mat' ...
                 'S1020V11_AFC_RightACL0_2409171645.mat' ...                  
                 };            
elseif subj==16
    filenames = {
                 'S1026V9_AFC_RightACL0_2409251353.mat' ...
                 'S1026V9_AFC_RightACL0_2409251346.mat' ...
                 'S1026V9_AFC_RightACL0_2409251339.mat' ...
                 'S1026V9_AFC_RightACL0_2409251331.mat' ...
                 'S1026V9_AFC_RightACL0_2409251322.mat' ... 
                 'S1026V9_AFC_RightACL0_2409251315.mat' ...   
                 'S1026V9_AFC_RightACL0_2409251308.mat' ... 
                 'S1026V9_AFC_RightACL0_2409251300.mat' ...     
                 'S1026V9_AFC_RightACL0_2409251253.mat' ...     
                 'S1026V9_AFC_RightACL0_2409251246.mat' ...    
                 'S1026V9_AFC_RightACL0_2409261655.mat' ...     
                 'S1026V9_AFC_RightACL0_2409261649.mat' ...     
                 'S1026V9_AFC_RightACL0_2409261643.mat' ...                     
                 }; 
elseif subj==5
    filenames = {
                 'S1015V10_AFC_RightACL0_2409261217.mat' ...
                 'S1015V10_AFC_RightACL0_2409261208.mat' ...
                 'S1015V10_AFC_RightACL0_2409261202.mat' ...
                 'S1015V10_AFC_RightACL0_2409261155.mat' ...
                 'S1015V10_AFC_RightACL0_2409261149.mat' ... 
                 'S1015V10_AFC_RightACL0_2409261143.mat' ...   
                 'S1015V10_AFC_RightACL0_2409261136.mat' ... 
                 'S1015V10_AFC_RightACL0_2409261130.mat' ...     
                 'S1015V10_AFC_RightACL0_2409261122.mat' ...     
                 'S1015V10_AFC_RightACL0_2410031542.mat' ... 
                 'S1015V10_AFC_RightACL0_2410031551.mat' ...     
                 'S1015V10_AFC_RightACL0_2410031559.mat' ...                      
                 };     
elseif subj==17
    filenames = {
                 'S1027V9_AFC_RightACL0_2410031133.mat' ...
                 'S1027V9_AFC_RightACL0_2410031126.mat' ...
                 'S1027V9_AFC_RightACL0_2410031110.mat' ...
                 'S1027V9_AFC_RightACL0_2410031102.mat' ...
                 'S1027V9_AFC_RightACL0_2410031049.mat' ... 
                 'S1027V9_AFC_RightACL0_2410031041.mat' ...   
                 'S1027V9_AFC_RightACL0_2410031033.mat' ... 
                 'S1027V9_AFC_RightACL0_2410031023.mat' ... 
                 'S1027V9_AFC_RightACL0_2410081020.mat' ... 
                 'S1027V9_AFC_RightACL0_2410081027.mat' ... 
                 'S1027V9_AFC_RightACL0_2410081033.mat' ... 
                 'S1027V9_AFC_RightACL0_2410081045.mat' ... 
                 };             
elseif subj==18
    filenames = {
                 'S1028V9_AFC_RightACL0_2411011027.mat' ...
                 'S1028V9_AFC_RightACL0_2411011020.mat' ...
                 'S1028V9_AFC_RightACL0_2411011011.mat' ...
                 'S1028V9_AFC_RightACL0_2411011035.mat' ...
                 'S1028V9_AFC_RightACL0_2411011043.mat' ...
                 'S1028V9_AFC_RightACL0_2411011050.mat' ...
                 'S1028V9_AFC_RightACL0_2411011107.mat' ...
                 'S1028V9_AFC_RightACL0_2411011114.mat' ...
                 'S1028V9_AFC_RightACL0_2411051339.mat' ...
                 'S1028V9_AFC_RightACL0_2411051327.mat' ...
                 'S1028V9_AFC_RightACL0_2411051313.mat' ...
                 };                 
elseif subj==20
    filenames = {
                 'S1030V8_AFC_RightACL0_2411041452.mat' ...
                 'S1030V8_AFC_RightACL0_2411041441.mat' ...
                 'S1030V8_AFC_RightACL0_2411041426.mat' ...
                 'S1030V8_AFC_RightACL0_2411041416.mat' ...
                 'S1030V8_AFC_RightACL0_2411151008.mat' ...
                 'S1030V8_AFC_RightACL0_2411151021.mat' ...
                 'S1030V8_AFC_RightACL0_2411151031.mat' ...
                 'S1030V8_AFC_RightACL0_2411151041.mat' ...
                 'S1030V8_AFC_RightACL0_2411151052.mat' ...
                 'S1030V8_AFC_RightACL0_2411151102.mat' ...
                 'S1030V8_AFC_RightACL0_2411151112.mat' ...
                 };                     
end

rgb = []; % COLOR CONDITIONS OF ACCOMMODATIVE STIMULUS--ALWAYS THE SAME
meanFocstmOptDst = []; % ACCOMMODATIVE STIMULUS DISTANCE
focStmOptDstIncr = []; % ACUITY STIMULUS DISTANCE
rspAcu = []; % RESPONSES
stimOrientation = []; % STIM ORIENTATION: 1 FOR -15 DEG, 2 FOR 15 DEG
indAcuRB = []; % 1 FOR RED, 2 FOR GREEN, 3 FOR BLUE

% FINDING PEAKS (AND PLOTTING IF NECESSARY)
for i = 1:length(filenames)
    load(fullfile(dataDirectory,filenames{i})); % LOAD FILE
    % STACK UP VARIABLES REQUIRED FOR ANALYSIS INTO COLUMNS
    rgb = [rgb; AFCp.rgb]; 
    meanFocstmOptDst = [meanFocstmOptDst; AFCp.meanFocstmOptDst];
    focStmOptDstIncr = [focStmOptDstIncr; AFCp.focStmOptDstIncr];
    rspAcu = [rspAcu; AFCp.rspAcu'];
    stimOrientation = [stimOrientation; AFCp.stimOrientation];
    indAcuRB = [indAcuRB; AFCp.indAcuRGBall(1:end-1)];
end

% COLORS FOR PLOTTING
rgbUnq = [1 0 0; 0 1 0; 0 0 1];

scaleFac = 0.816; % ACCOUNT FOR OPTICAL PROPERTIES OF BVAMS

if bPLOT
    figure;
end

defocusLCAmeasured = []; % BEST DISTANCE
dprimeFitAll = {}; % D PRIMES
for rgbAcuCnd = 1:3 % SORT rgbAcuCnd DATA BY COLOR
    unqFocDst = unique(focStmOptDstIncr(indAcuRB==rgbAcuCnd));
    for i = 1:length(unqFocDst) % LOOP OVER STIMULUS DISTANCE
        indAnalysis = abs(focStmOptDstIncr-unqFocDst(i))<0.001 & indAcuRB==rgbAcuCnd;
        % PERCENT CORRECT
        PC(i) = sum(rspAcu(indAnalysis)==stimOrientation(indAnalysis))./sum(indAnalysis);
        % 95% CONFIDENCE INTERVALS
        PCci(:,i) = binoinv([0.025 0.975],sum(focStmOptDstIncr==unqFocDst(i) & indAcuRB==rgbAcuCnd),PC(i))./sum(focStmOptDstIncr==unqFocDst(i) & indAcuRB==rgbAcuCnd);
    end

    % FITTING SPLINE TO DATA
    PCfitSupport = min(unqFocDst.*scaleFac):0.01:max(unqFocDst.*scaleFac);
    PCfit = spline(unqFocDst.*scaleFac,PC,PCfitSupport);
    if rgbAcuCnd==2 || rgbAcuCnd==3 % IF GREEN OR BLUE CONDITIONS
        % PROTECT AGAINST IMPLAUSIBLE PEAK VALUES FROM MEASUREMENT NOISE OR
        % MONOCHROMATIC ABERRATIONS BY ENFORCING CONSTRAINT THAT PEAK FOR 
        % GREEN IS WITHIN 1.0D OF PEAK FOR RED, AND PEAK FOR BLUE IS WITHIN 
        % 1.0D OF PEAK FOR GREEN. THESE LIMITS ARE MUCH LARGER THAN ANY LCA
        % MAGNITUDE OBSERVED IN ANY STUDY. ONLY MAKES A DIFFERENCE FOR
        % SUBJECT S5. 
        indValid = abs(PCfitSupport-defocusLCAmeasured(rgbAcuCnd-1))<1;
        PCfitSupportValid = PCfitSupport(indValid);
        PCfitValid = PCfit(indValid);
        % FIND PEAK
        [~,indLCA] = max(PCfitValid);
        defocusLCAmeasured(rgbAcuCnd) = PCfitSupportValid(indLCA);   
    else
        % FIND PEAK
        [~,indLCA] = max(PCfit);
        defocusLCAmeasured(rgbAcuCnd) = PCfitSupport(indLCA);        
    end
    if bPLOT
        hold on;
        plot(PCfitSupport,PCfit,'-','Color',rgbUnq(rgbAcuCnd,:),'MarkerFaceColor','w','LineWidth',2,'MarkerSize',20);
        plot(unqFocDst.*scaleFac,PC,'o','Color',rgbUnq(rgbAcuCnd,:),'MarkerFaceColor','w','LineWidth',2,'MarkerSize',20);
        errorbar(unqFocDst.*scaleFac,PC,PC-PCci(1,:),PCci(2,:)-PC,'o','Color',rgbUnq(rgbAcuCnd,:),'MarkerFaceColor','w','LineWidth',2,'MarkerSize',20);
        ylim([0.4 1]);
        if rgbAcuCnd==1
           formatFigure('Relative optical distance (D)','Proportion Correct',['Subject ' num2str(subj)]);
        else
           formatFigure('Relative optical distance (D)','Proportion Correct'); 
        end
    end
    % SUBTRACT EPSILON IF SUBJECT HITS PERFORMANCE OF 100%
    epsilonPC = 0.99;
    PCfit(PCfit>epsilonPC) = epsilonPC;
    % CALCULATE D-PRIME
    dprimeFitAll{rgbAcuCnd} = 2*norminv(PCfit);
end

% THE BLOCK OF CODE BELOW DOES BOOTSTRAPPING--NEARLY IDENTICAL TO CODE
% ABOVE, JUST WITH MANY LOOPS FOR BOOTSTRAPS
defocusLCAmeasuredBoots = []; % BOOTSTRAPPING BEST DISTANCE VALUES
if nBoots>1
    % INITIALIZE TEMPORARY ARRAY FOR STORING BEST DISTANCE VALUES
    defocusLCAmeasuredTmp = [];
    for j = 1:nBoots
        for rgbAcuCnd = 1:3 % SORT rgbAcuCnd DATA BY COLOR
            unqFocDstBoots = unique(focStmOptDstIncr(indAcuRB==rgbAcuCnd));
            for i = 1:length(unqFocDstBoots) % LOOP OVER STIMULUS DISTANCE
                indAnalysisBoots = abs(focStmOptDstIncr-unqFocDstBoots(i))<0.001 & indAcuRB==rgbAcuCnd;
                % COMPUTE PERCENT CORRECT
                indBoots = randsample(find(indAnalysisBoots),sum(indAnalysisBoots),'true');                
                % PERCENT CORRECT
                PCboots(i) = sum(rspAcu(indBoots)==stimOrientation(indBoots))./sum(indBoots);
            end
        
            % FITTING SPLINE TO DATA
            PCfitSupportBoots = min(unqFocDstBoots.*scaleFac):0.01:max(unqFocDstBoots.*scaleFac);
            PCfitBoots = spline(unqFocDstBoots.*scaleFac,PCboots,PCfitSupportBoots);
            if rgbAcuCnd==2 || rgbAcuCnd==3 % IF GREEN OR BLUE CONDITIONS
                % PROTECT AGAINST IMPLAUSIBLE PEAK VALUES FROM MEASUREMENT NOISE OR
                % MONOCHROMATIC ABERRATIONS BY ENFORCING CONSTRAINT THAT PEAK FOR 
                % GREEN IS WITHIN 1.0D OF PEAK FOR RED, AND PEAK FOR BLUE IS WITHIN 
                % 1.0D OF PEAK FOR GREEN. THESE LIMITS ARE MUCH LARGER THAN ANY LCA
                % MAGNITUDE OBSERVED IN ANY STUDY. ONLY MAKES A DIFFERENCE FOR
                % SUBJECT S5.
                indValidBoots = abs(PCfitSupportBoots-defocusLCAmeasuredTmp(rgbAcuCnd-1))<1;
                PCfitSupportValidBoots = PCfitSupportBoots(indValidBoots);
                PCfitValidBoots = PCfitBoots(indValidBoots);
                % FIND PEAK
                [~,indLCAboots] = max(PCfitValidBoots);
                defocusLCAmeasuredTmp(rgbAcuCnd) = PCfitSupportValidBoots(indLCAboots);   
            else
                % FIND PEAK
                [~,indLCAboots] = max(PCfitBoots);
                defocusLCAmeasuredTmp(rgbAcuCnd) = PCfitSupportBoots(indLCAboots);        
            end
        end
        % STORE BOOTSTRAPPED PEAK DEFOCUS VALUES
        defocusLCAmeasuredBoots(:,j) = defocusLCAmeasuredTmp;
    end
end

wavePlotLCA = 380:875; % SUPPORT FOR LCA FUNCTION
wavePlotPrimaries = [616 533 468]; % PEAK WAVELENGTHS OF PRIMARIES
nRepeatFit = 100; % NUMBER OF TIMES TO REPEAT FIT
% LCA PARAMETERS
q1all = [];
q2all = [];
q3all = [];
errAll = [];
% REPEATING FIT OF LCA FUNCTION TO ACCOUNT FOR INSTABILITIES
for i = 1:nRepeatFit
   [q1,q2,q3,err] = ARC_LCAfit(wavePlotPrimaries,-defocusLCAmeasured);
   q1all(i) = q1;
   q2all(i) = q2;
   q3all(i) = q3;
   errAll(i) = err;
   if mod(i,1000)==0
       display(['Iteration ' num2str(i)]);
   end
end
% GET BEST OF REPEATED FITS
[~,indBestFit] = min(errAll);
q1best = q1all(indBestFit);
q2best = q2all(indBestFit);
q3best = q3all(indBestFit);
% DEFOCUS AT GREEN (NOT CURRENTLY USING BUT COULD FOR OTHER PURPOSES)
Dgreen = humanWaveDefocusParameterized(532,q1best,q2best,q3best);

if bPLOT
    figure; 
    hold on;
    plot(wavePlotLCA,humanWaveDefocusParameterized(wavePlotLCA,q1best,q2best,q3best),'k-','LineWidth',1); 
    plot(wavePlotPrimaries,-defocusLCAmeasured,'ko','MarkerSize',15,'MarkerFaceColor','w');
    axis square; 
    set(gca,'FontSize',15); 
    xlabel('Wavelength (\lambda)'); 
    ylabel('Relative Defocus (D)');
    xlim([400 875]);
end

end
