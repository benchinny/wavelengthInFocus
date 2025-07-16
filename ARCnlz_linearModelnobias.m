function [weightsRBS1_x, weightsRBS1_y, rhoFull, rhoNoColor, rhoColor, aic, aicNoColor, weightsRBSci, trialMeans, meanChanges, deltaR, deltaB, deltaS] = ARCnlz_linearModelnobias(sn,bPLOT,nBoot,bLOWLUM)

bEXCLUDE = true;
gammaFactorR = 2.4;
gammaFactorB = 2.2;
scaleFactor = 0.8;
bRELATIVE_LUM = 0;
maxLumCdm2 = 0.40;

if sn==11 % 'VISIT' NUMBERS
   % vs = [2 3 4 5 6 7];
   vs = [5 6];
   excludeTrials = [];
elseif sn==12
   vs = [3:7];
   excludeTrials = [];
elseif sn==13
   vs = 1:4;
   excludeTrials = [];
elseif sn==14
   vs = 1:4;
   excludeTrials = [];    
elseif sn==15
   vs = 7:12;
   excludeTrials = [];  
elseif sn==16
   % vs = 7:12;
   vs = [11 12];  
   excludeTrials = [];     
elseif sn==17
   % vs = 1:6;
   vs = [5 6];
   excludeTrials = [];  
elseif sn==18
   vs = 7:10;
   % vs = [5 6];
   excludeTrials = [];     
elseif sn==19
   % vs = 5:10;
   vs = [9 10];
   excludeTrials = [];     
elseif sn==21
   % vs = 1:6;
   vs = [5 6];
   excludeTrials = [];
elseif sn==23
   % vs = 14:21;
   vs = 18:21;
   excludeTrials = [];   
elseif sn==24
   % vs = 9:16;
   vs = 13:16;
   excludeTrials = [];     
elseif sn==25
   % vs = 6:13;
   vs = 10:13;
   excludeTrials = [];        
elseif sn==26
   % vs = 4:11;
   vs = 8:11;
   excludeTrials = [];  
elseif sn==27
   % vs = 4:11;
   vs = 9:12;
   excludeTrials = [];     
elseif sn==28
   % vs = 4:11;
   vs = 17:20;
   excludeTrials = [];              
elseif sn==29
   % vs = 4:11;
   vs = 11:14;
   excludeTrials = [];        
elseif sn==30
   % vs = 4:11;
   vs = 13:16;
   excludeTrials = [];       
elseif sn==31
   % vs = 4:11;
   vs = 11:14;
   excludeTrials = [];          
elseif sn==32
   % vs = 4:11;
   vs = 15:18;
   excludeTrials = [];           
elseif sn==33
   % vs = 4:11;
   vs = 8:11;
   excludeTrials = []; 
elseif sn==7
   vs = 10:14;
   excludeTrials = [];
else
   error('ARCnlz_linearModelnobias: unhandled subject number!');
end
ey = 1; % 1 for Right eye, 2 for Binocular ');
filePath = 'G:\My Drive\exp_bvams\code_repo\';


if strcmp(getenv('username'),'bankslab')
   dataDirectory = [filePath 'ARC\'];
elseif strcmp(getenv("USER"),'benchin')
   dataDirectory = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Analysis/'; 
elseif strcmp(getenv("USER"),'emily')
   dataDirectory = '/Users/emily/Library/CloudStorage/GoogleDrive-emilyacooper@gmail.com/Shared drives/ARChroma/Analysis/';
elseif strcmp(getenv("USER"),'benjaminchin')
   dataDirectory = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Analysis/';
elseif strcmp(getenv("USER"),'ben')
   dataDirectory = '/home/ben/Documents/ARchroma/Analysis/';   
end

% THIS JUST LOADS A FILE CONTAINING FILE NAMES OF .mat FILES CONTAINING
% METADATA FOR EACH EXPERIMENT BLOCK
if ey(1)==2
    load([dataDirectory 'AFCflsB.mat']); c1=2;
elseif ey(1)==1
    load([dataDirectory 'AFCflsR.mat']); c1=1;
end
x = [];
y = [];
t3 = [];
for i = 1:length(vs)
    if strcmp(getenv('username'),'bankslab')
        fnm=AFCfls{sn, vs(i)};
        fnmTmp = fnm;
        load(fnm);
    else
        % LOADS .mat FILE CONTAINING METADATA FOR EXPERIMENT BLOCK
        fnmTmp=AFCfls{sn, vs(i)};         
        load([dataDirectory fnmTmp(37:end)]);
    end
    if i==1
       AFCp = rmfield(AFCp,'v0');
       AFCp = rmfield(AFCp,'AFCv');
       AFCpAll = AFCp;
    else
       AFCp = rmfield(AFCp,'v0'); 
       AFCp = rmfield(AFCp,'AFCv');
       AFCpAll = structmerge(AFCpAll,AFCp,length(AFCp.v00));
    end
    dateCodeAll = [];
    % CONVERTING DATE OF .mat FILE TO JSON FILE NAME
    dateCode = fnmTmp((end-9):end);
    dateCodeAll = [dateCodeAll; dateCode];
    dateCode(end) = dateCode(end)-1; % THE NEXT FEW LINES ARE FOR ADDING ROBUSTNESS
    dateCodeAll = [dateCodeAll; dateCode];
    dateCode(end) = dateCode(end)+2;
    dateCodeAll = [dateCodeAll; dateCode];
    if strcmp(getenv('username'),'bankslab')
        jsonDirectory = 'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real\';
    else
        jsonDirectory = dataDirectory;
    end    
    jsonFile = [];
    for j = 1:size(dateCodeAll,1) % FOR ROBUSTNESS: SOMETIMES THE FILENAMES DON'T MATCH EXACTLY
        jsonFileStr = ['20' dateCodeAll(j,1:2) '-' dateCodeAll(j,3:4) '-' dateCodeAll(j,5:6) ' ' dateCodeAll(j,7:8) '.' dateCodeAll(j,9:10) '.json'];
        if isfile([jsonDirectory jsonFileStr])
            jsonFile = jsonFileStr;
        end
    end
    jsonPath = jsonDirectory;
    dt=jsondecode(fileread([jsonPath jsonFile]));
    % GRABS RAW PIXEL SEPARATIONS BETWEEN AUTOREFRACTOR BARS
    x=[x; -1.*[(cell2mat(struct2cell(dt.ext_right_mu))-cell2mat(struct2cell(dt.ext_left_mu)))]];
    y=[y; -1.*[(cell2mat(struct2cell(dt.ext_bottom_mu))-cell2mat(struct2cell(dt.ext_top_mu)))]];

    % THIS BLOCK SORTS THE DATA BY TRIAL
    t3=[t3; cell2mat(struct2cell(dt.time_totalsecs))];           
end
AFCp = AFCpAll;

startIndices = [1; find(diff(t3)>1)+1];
startTimesFromAutoref = t3(startIndices);
startTimesFromComputer = squeeze(AFCp.t3(:,6,1))+squeeze(AFCp.t3(:,5,1))*60+squeeze(AFCp.t3(:,4,1))*3600;
endIndices = [startIndices(2:end)-1; length(t3)];
endTimesFromAutoref = t3(endIndices);
endTimesFromComputer = squeeze(AFCp.t3(:,6,2))+squeeze(AFCp.t3(:,5,2))*60+squeeze(AFCp.t3(:,4,2))*3600;

v0=[find(diff(t3)>1); length(t3)];
i0=1;
for k0=1:length(v0);
   v1=i0:v0(k0);
   x0{k0}=x(v1);
   y0{k0}=y(v1);
   i0=v0(k0)+1;
end

% FINDS AUTOREFRACTOR VALUES FOR 'REFERENCE' ACCOMMODATION (SOME
% FLEXIBILITY HERE ABOUT HOW TO DEFINE THE REFERENCE)
i0=find(AFCp.v00(:,1)==0);
i0 = 1;
x1=[]; y1=[];
for k0=1:length(i0)
 x1=[x1 x0{i0(k0)}'];
 y1=[y1 y0{i0(k0)}'];  
end
x1 = x1(1:10);
y1 = y1(1:10);

% MAIN ANALYSIS
uniqueConditions = unique([AFCp.rgb100 AFCp.rgb200 AFCp.v00],'rows');
uniqueRGBvalues = unique([AFCp.rgb100 AFCp.rgb200],'rows');
xScale = -2.87; % SCALE FACTOR CONVERTING PIXELS TO DIOPTERS
yScale = -2.58; % SCALE FACTOR CONVERTING PIXELS TO DIOPTERS
optDistScale = 0.8;
timeSeries = {};
rgbValuesAll = {};
x3stack = zeros([length(AFCp.v00) 400]);
y3stack = zeros([length(AFCp.v00) 400]);
meanChangeXvec = [];
meanChangeYvec = [];

for i = 1:size(uniqueConditions,1)
    % INDEX OF COLOR CONDITIONS
    indCnd = find(ismember([AFCp.rgb100 AFCp.rgb200 AFCp.v00],uniqueConditions(i,:),'rows'));
    if bEXCLUDE
       indCnd = indCnd(~ismember(indCnd,excludeTrials));
    end
    clear x2 y2 x3 y3;
    x3=[]; y3=[]; 
    rgbValues = [];
    trialMarkerForPlot = [];
    meanChangeXtmp = [];
    meanChangeYtmp = [];
    for k0=1:size(indCnd,1)
       % ANALYZING SUBJECT'S ACCOMMODATION
       x2{k0}=x0{indCnd(k0)}-mean(x1); % MEAN CENTERING
       y2{k0}=y0{indCnd(k0)}-mean(y1); % MEAN CENTERING
       x3tmp = (x2{k0}')./xScale;
       y3tmp = (y2{k0}')./yScale;
       % REMOVING OUTLIERS
       x3diff = [0 diff(x3tmp)];
       y3diff = [0 diff(y3tmp)];
       x3outliers = abs(x3diff)>1 | abs(x3tmp)>5;
       y3outliers = abs(y3diff)>1 | abs(y3tmp)>5;
       meanx3 = mean(x3tmp(~x3outliers));
       meany3 = mean(y3tmp(~y3outliers));
       if bEXCLUDE
           x3tmp(x3outliers) = meanx3;
           y3tmp(y3outliers) = meany3;
       end       
       x3=[x3 x3tmp];
       y3=[y3 y3tmp];
       x3stack(indCnd(k0),1:length(x3tmp)) = x3tmp;
       y3stack(indCnd(k0),1:length(y3tmp)) = y3tmp;
       trialMarkerForPlot(k0) = length(x3);
       % ACCOMMODATIVE DEMAND FROM EXPERIMENT 
       sinValuesTmp = AFCp.sinValues(indCnd(k0),:);
       tSin = 0:(1/(length(sinValuesTmp)-1)):1;
       tSinInterp = 0:(1/(length(x2{k0})-1)):1;
       accContinuous = interp1(tSin,sinValuesTmp,tSinInterp); 
       % THIS IS AN OBNOXIOUS WAY OF COMPUTING THE AVERAGE CHANGE WITHIN A
       % TRIAL, BUT IT WORKS AND MAY BE A BIT MORE ROBUST
%       diffVec = imresize([-2/length(accContinuous) 2/length(accContinuous)],size(accContinuous),'nearest');
       diffVec = imresize([0 -4/length(accContinuous) 0 4/length(accContinuous)],size(accContinuous),'nearest');
%        if abs(corr(accContinuous',diffVec'))<0.95
%            error('ARCnlz_linearModelnobias: you may want to check whether the step change occurs halfway through the trial, or not!');
%        end
%        meanChangeX{i} = sum(diffVec.*x3tmp.*xScale)./xScale;
%        meanChangeY{i} = sum(diffVec.*y3tmp.*yScale)./yScale;
       meanChangeXtmp(k0) = sum(diffVec.*x3tmp.*xScale)./xScale;
       meanChangeYtmp(k0) = sum(diffVec.*y3tmp.*yScale)./yScale;
       % VECTOR OF RGB VALUES FOR PLOTTING
       rgbValues = [rgbValues imresize([AFCp.rgb100(indCnd(k0),:)' AFCp.rgb200(indCnd(k0),:)'],[3 length(tSinInterp)],'nearest')];
    end

    % STORE ACCOMMODATION AND COLOR VALUES FOR PLOTTING
    timeSeries{i,1} = x3;
    timeSeries{i,2} = y3;
    rgbValuesAll{i} = rgbValues;
    trialMarkerForPlotCell{i} = trialMarkerForPlot;
    indCndCell{i} = indCnd;
    meanChangeXvec(indCnd) = meanChangeXtmp;
    meanChangeYvec(indCnd) = meanChangeYtmp;
end

scaleEquateRB = 4;

if bRELATIVE_LUM
   deltaR = scaleEquateRB.*(AFCp.rgb200(:,1).^gammaFactorR)./(scaleEquateRB.*(AFCp.rgb200(:,1).^gammaFactorR) + (AFCp.rgb200(:,3).^gammaFactorB)) ...
             - scaleEquateRB.*(AFCp.rgb100(:,1).^gammaFactorR)./(scaleEquateRB.*(AFCp.rgb100(:,1).^gammaFactorR) + (AFCp.rgb100(:,3).^gammaFactorB));
   deltaB = (AFCp.rgb200(:,3).^gammaFactorB)./((AFCp.rgb200(:,3).^gammaFactorB) + scaleEquateRB.*(AFCp.rgb200(:,1).^gammaFactorR)) ...
             - (AFCp.rgb100(:,3).^gammaFactorB)./((AFCp.rgb100(:,3).^gammaFactorB) + scaleEquateRB.*(AFCp.rgb100(:,1).^gammaFactorR));
   deltaB = [];
else
   deltaR = scaleEquateRB.*AFCp.rgb200(:,1).^gammaFactorR - scaleEquateRB.*AFCp.rgb100(:,1).^gammaFactorR;
   deltaB = AFCp.rgb200(:,3).^gammaFactorB - AFCp.rgb100(:,3).^gammaFactorB;
end
% COMPONENTS OF LINEAR REGRESSION
deltaS = AFCp.v00*scaleFactor;
delta1 = ones(size(deltaR));
deltaR = deltaR.*maxLumCdm2;
deltaB = deltaB.*maxLumCdm2;

% IF EXCLUDING HIGH LUMINANCE TRIALS
lumCutoff = 2;
if bLOWLUM
    indLowLum = 0.4.*(scaleEquateRB.*AFCp.rgb100(:,1).^gammaFactorR + AFCp.rgb100(:,3).^gammaFactorB)<lumCutoff | ...
                0.4.*(scaleEquateRB.*AFCp.rgb200(:,1).^gammaFactorR + AFCp.rgb200(:,3).^gammaFactorB)<lumCutoff;
    % indLowLum = ~indLowLum;
    deltaR = deltaR(indLowLum);
    deltaB = deltaB(indLowLum);
    deltaS = deltaS(indLowLum);
    meanChangeXvec = meanChangeXvec(indLowLum);
    meanChangeYvec = meanChangeYvec(indLowLum);
end
meanChanges = meanChangeXvec;

% DOING THE LINEAR REGRESSION
weightsRBS1_x = [deltaR deltaB deltaS]\(meanChangeXvec');
weightsRBS1_y = [deltaR deltaB deltaS]\(meanChangeYvec');
weightsS1_x = [deltaS]\(meanChangeXvec');
weightsS1_y = [deltaS]\(meanChangeYvec');
weightsRB_x = [deltaR deltaB]\(meanChangeXvec');
weightsRB_y = [deltaR deltaB]\(meanChangeYvec');

% BOOTSTRAPPING
for i = 1:nBoot
   indBoot = randsample(1:length(meanChangeXvec),length(meanChangeXvec),true);
   if bRELATIVE_LUM
       weightsRBSboot(i,:) = [deltaR(indBoot) deltaS(indBoot)]\(meanChangeXvec(indBoot)');
   else
       weightsRBSboot(i,:) = [deltaR(indBoot) deltaB(indBoot) deltaS(indBoot)]\(meanChangeXvec(indBoot)');
   end
end
weightsRBSci = quantile(weightsRBSboot,[0.16 0.84]);

% COMPUTING CORRELATIONS BETWEEN DIFFERENT MODEL PREDICTIONS AND DATA
rhoFull = corr([deltaR deltaB deltaS]*weightsRBS1_x,meanChangeXvec');
rhoNoColor = corr([deltaS]*weightsS1_x,meanChangeXvec');
rhoColor = corr([deltaR deltaB]*weightsRB_x,meanChangeXvec');

% COMPUTING CORRELATIONS BETWEEN DIFFERENT MODEL PREDICTIONS AND DATA
% rhoFull = corr([deltaR deltaB deltaS]*weightsRBS1_y,meanChangeYvec');
% rhoNoColor = corr([deltaS]*weightsS1_y,meanChangeYvec');
% rhoColor = corr([deltaR deltaB]*weightsRB_y,meanChangeYvec');

if bRELATIVE_LUM
    nParams = 3;
    nParamsNoColor = 2;
else
    nParams = 4;
    nParamsNoColor = 2;
end

% COMPUTING COMPONENETS FOR AIC
trialMeans = [deltaR deltaB deltaS]*weightsRBS1_x;
errorIndividual = meanChangeXvec' - trialMeans;
for i = 1:100
   [stdTmp(i),LLtmp(i)] = ARCfitStdGauss(errorIndividual);
end
[~,bestInd] = min(LLtmp);
estResidualStd = stdTmp(bestInd);
LL = sum(log(normpdf(meanChangeXvec',trialMeans,estResidualStd)));
aic = 2*nParams-2*LL;

% COMPUTING COMPONENTS FOR AIC FOR NO-COLOR MODEL
trialMeansNoColor = [deltaS]*weightsS1_x;
errorIndividualNoColor = meanChangeXvec' - trialMeansNoColor;
for i = 1:100
   [stdTmp(i),LLtmp(i)] = ARCfitStdGauss(errorIndividualNoColor);
end
[~,bestInd] = min(LLtmp);
estResidualStdNoColor = stdTmp(bestInd);
LLnoColor = sum(log(normpdf(meanChangeXvec',trialMeansNoColor,estResidualStdNoColor)));
aicNoColor = 2*nParamsNoColor-2*LLnoColor;

% FOR COMPARING MANUAL AIC COMPUTATION WITH BUILT-IN MATLAB FUNCTION
if bRELATIVE_LUM
    % put data into table in prep for running built in matlab models
    tbl = table(deltaR,deltaS,meanChangeXvec','VariableNames',{'deltaR','deltaS','resp'});
    
    % fit two linear mixture models to data, one with color and one
    % without, run model comparison
    lme_nocolor = fitlme(tbl,'resp~-1+deltaS');
    lme = fitlme(tbl,'resp~-1+deltaS+deltaR');
    compare(lme,lme_nocolor);
    
%    aic = 2*lme.NumPredictors - 2*lme.LogLikelihood;
%    aicNoColor = 2*lme_nocolor.NumPredictors - 2*lme_nocolor.LogLikelihood;
else
    % put data into table in prep for running built in matlab models
    tbl = table(deltaR,deltaB,deltaS,meanChangeXvec','VariableNames',{'deltaR','deltaB','deltaS','resp'});
    
    % fit two linear mixture models to data, one with color and one
    % without, run model comparison
    lme_nocolor = fitlme(tbl,'resp~-1+deltaS');
    lme = fitlme(tbl,'resp~-1+deltaS+deltaB+deltaR');
    compare(lme,lme_nocolor);
    
%    aic = 2*lme.NumPredictors - 2*lme.LogLikelihood;
%    aicNoColor = 2*lme_nocolor.NumPredictors - 2*lme_nocolor.LogLikelihood;
end

if bPLOT
    figure;
    set(gcf,'Position',[189 395 1280 420]);
    subplot(1,3,1);
%    plot([deltaR deltaB deltaS]*weightsRBS1_x,meanChangeXvec,'ko','LineWidth',1);
    predictionTmp = [deltaR deltaB deltaS]*weightsRBS1_x;
    hold on;
    for i = 1:length(predictionTmp)
        deltaRtmp = (1/maxLumCdm2)*deltaR(i);
        deltaBtmp = (1/maxLumCdm2)*deltaB(i);
        RBratio = 0.25.*(deltaRtmp-deltaBtmp+2);
        RBratio
        plot(predictionTmp(i),meanChangeXvec(i),'o','LineWidth',1,'Color',[RBratio 0 1-RBratio],'MarkerFaceColor',[RBratio 0 1-RBratio]);
        % plot(predictionTmp(i),meanChangeXvec(i),'o','LineWidth',1,'Color','k','MarkerFaceColor','w');
    end    
    plot([-2 2],[-2 2],'k--');
    xlim([-2 2]);
    ylim([-2 2]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaD');
    ylabel('Measured \DeltaD');
    title(['Correlation = ' num2str(rhoFull,3)]);
    axis square;
    
    subplot(1,3,2);
    plot([deltaR deltaB deltaS]*weightsRBS1_y,meanChangeYvec,'ko','LineWidth',1);
    xlim([-2 2]);
    ylim([-2 2]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaD');
    ylabel('Measured \DeltaD');
    axis square;
    
    subplot(1,3,3);
    hold on;
    bar(1,weightsRBS1_x(1),'FaceColor','r');
    if ~bRELATIVE_LUM
       bar(2,weightsRBS1_x(2),'FaceColor','b');
       bar(3,weightsRBS1_x(3),'FaceColor','k');
       set(gca,'XTick',[1 2 3]);
       set(gca,'XTickLabel',{'Red' 'Blue' 'D_{opt}'});
    else
       bar(2,weightsRBS1_x(2),'FaceColor','k');
       set(gca,'XTick',[1 2]);
       set(gca,'XTickLabel',{'R-B ratio' 'D_{opt}'}); 
    end
    title('Weights');
    set(gca,'FontSize',20);
    ylim(max(weightsRBS1_x).*[-1.2 1.2]);
    axis square;
    
    figure;
    set(gcf,'Position',[189 395 1280 420]);
    subplot(1,3,1);
    % plot([deltaS]*weightsS1_x,meanChangeXvec,'ko','LineWidth',1);
    predictionTmp = [deltaS]*weightsS1_x;
    hold on;
    for i = 1:length(predictionTmp)
        deltaRtmp = (1/maxLumCdm2)*deltaR(i);
        deltaBtmp = (1/maxLumCdm2)*deltaB(i);
        RBratio = 0.25.*(deltaRtmp-deltaBtmp+2);
        RBratio
        plot(predictionTmp(i),meanChangeXvec(i),'o','LineWidth',1,'Color',[RBratio 0 1-RBratio],'MarkerFaceColor',[RBratio 0 1-RBratio]);
        % plot(predictionTmp(i),meanChangeXvec(i),'o','LineWidth',1,'Color','k','MarkerFaceColor','w');
    end
    plot([-2 2],[-2 2],'k--');
    xlim([-2 2]);
    ylim([-2 2]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaA');
    ylabel('Measured \DeltaA');
    title(['Correlation = ' num2str(rhoNoColor,3)]);
    axis square;
    
    subplot(1,3,2);
    plot([deltaS]*weightsS1_y,meanChangeYvec,'ko','LineWidth',1);
    xlim([-2 2]);
    ylim([-2 2]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaD');
    ylabel('Measured \DeltaD');
    axis square;
    
    subplot(1,3,3);
    hold on;
    bar(1,weightsS1_x(1),'FaceColor','k');
    set(gca,'XTick',[1 2]);
    set(gca,'XTickLabel',{'D_{opt}'});
    title('Weights');
    set(gca,'FontSize',20);
    ylim(max(weightsS1_x).*[-1.2 1.2]);
    axis square;
    
    figure;
    set(gcf,'Position',[189 395 1280 420]);
    subplot(1,3,1);
    plot([deltaR deltaB]*weightsRB_x,meanChangeXvec,'ko','LineWidth',1);
    xlim([-2 2]);
    ylim([-2 2]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaD');
    ylabel('Measured \DeltaD');
    title(['Correlation = ' num2str(corr([deltaR deltaB]*weightsRB_x,meanChangeXvec'),3)]);
    axis square;
    
    subplot(1,3,2);
    plot([deltaR deltaB]*weightsRB_y,meanChangeYvec,'ko','LineWidth',1);
    xlim([-2 2]);
    ylim([-2 2]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaD');
    ylabel('Measured \DeltaD');
    axis square;
    
    subplot(1,3,3);
    hold on;
    if ~bRELATIVE_LUM    
        bar(1,weightsRB_x(1),'FaceColor','r');
        bar(2,weightsRB_x(2),'FaceColor','b');
        set(gca,'XTick',[1 2]);
        set(gca,'XTickLabel',{'R' 'B'});
    else
        bar(1,weightsRB_x(1),'FaceColor','r');
        set(gca,'XTick',[1]);
        set(gca,'XTickLabel',{'R'});
    end
    
    title('Weights');
    set(gca,'FontSize',20);
%    ylim(max(weightsRB_x).*[-1.2 1.2]);
    axis square;
end

end