function [d, wS, rbThreshold, rhoSwitch, rhoColorSwitch, aic, trialMeans,meanChanges, deltaRB1, deltaRB2, deltaS] = ARCnlz_linearSwitching(sn,bPLOT,dMin)

bEXCLUDE = true;
gammaFactorR = 2.4;
gammaFactorB = 2.2;
scaleFactor = 0.8;
maxLumCdm2 = 0.87;

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
   vs = 5:10;
   % vs = [9 10];
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
   vs = 6:13;
   % vs = 10:13;
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
elseif sn==32
   % vs = 4:11;
   vs = 15:18;
   excludeTrials = [];           
elseif sn==33
   % vs = 4:11;
   vs = 8:11;
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
       % diffVec = imresize([-2/length(accContinuous) 2/length(accContinuous)],size(accContinuous),'nearest');
       diffVec = imresize([0 -4/length(accContinuous) 0 4/length(accContinuous)],size(accContinuous),'nearest');
       % if abs(corr(accContinuous',diffVec'))<0.95
       %     error('ARCnlz_linearModelnobias: you may want to check whether the step change occurs halfway through the trial, or not!');
       % end
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
meanChanges = meanChangeXvec;

scaleEquateRB = 1/0.25;

% deltaR = scaleEquateRB.*AFCp.rgb200(:,1).^gammaFactorR - scaleEquateRB.*AFCp.rgb100(:,1).^gammaFactorR;
% deltaB = AFCp.rgb200(:,3).^gammaFactorB - AFCp.rgb100(:,3).^gammaFactorB;
deltaS = AFCp.v00*scaleFactor;

% COMPUTE DIFFERENCES BETWEEN RED AND BLUE LUMINANCES
deltaRB1 = scaleEquateRB.*AFCp.rgb100(:,1).^gammaFactorR - AFCp.rgb100(:,3).^gammaFactorB;
deltaRB2 = scaleEquateRB.*AFCp.rgb200(:,1).^gammaFactorR - AFCp.rgb200(:,3).^gammaFactorB;
deltaRB1 = deltaRB1.*maxLumCdm2;
deltaRB2 = deltaRB2.*maxLumCdm2;

for i = 1:100
   [dTmp(i), wStmp(i), rbThresholdTmp(i), mseTmp(i)] = ARCnlzSwitching(meanChangeXvec',deltaRB1,deltaRB2,deltaS,dMin);
   [~,indBest] = min(mseTmp);
   d = dTmp(indBest);
   wS = wStmp(indBest);
   rbThreshold = rbThresholdTmp(indBest);
end

dOrig = d;
wSorig = wS;
rbThresholdOrig = rbThreshold;

c = zeros(size(deltaS));
wS = wS.*ones(size(deltaS));

c(deltaRB1<rbThreshold & deltaRB2>rbThreshold)=d;
c(deltaRB1>rbThreshold & deltaRB2<rbThreshold)=-d;

deltaApredicted = wS.*deltaS + c;

rhoSwitch = corr(deltaApredicted,meanChangeXvec');
rhoColorSwitch = corr(c,meanChangeXvec');

trialMeans = deltaApredicted;
errorIndividual = meanChangeXvec' - trialMeans;
% COMPUTE ERROR OF RESIDUALS USING LOG-LIKELIHOOD METHOD
for i = 1:100
   [stdTmp(i),LLtmp(i)] = ARCfitStdGauss(errorIndividual);
end
[~,bestInd] = min(LLtmp);
estResidualStd = stdTmp(bestInd);
LL = sum(log(normpdf(meanChangeXvec',trialMeans,estResidualStd)));
nParams = 3;
aic = 2*nParams-2*LL;

if bPLOT
    figure;
    set(gcf,'Position',[189 395 1280 420]);
    subplot(1,3,1);
    hold on;
    for i = 1:length(deltaApredicted)
        if c(i)>0
            clrPlot = [0.6 0 0];
        elseif c(i)<0
            clrPlot = [0 0 1];
        else
            clrPlot = [0.6 0 1];
        end
        plot(deltaApredicted(i),meanChangeXvec(i),'o','LineWidth',1,'Color',clrPlot,'MarkerFaceColor',clrPlot);
    end
    xlim([-3 3]);
    ylim([-3 3]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaD');
    ylabel('Measured \DeltaD');
    title(['Correlation = ' num2str(rhoSwitch,3)]);
    axis square;
    % 
    % subplot(1,3,2);
    % plot([deltaR deltaB deltaS delta1]*weightsRBS1_y,meanChangeYvec,'ko','LineWidth',1);
    % xlim([-3 3]);
    % ylim([-3 3]);
    % set(gca,'FontSize',15);
    % xlabel('Prediction \DeltaD');
    % ylabel('Measured \DeltaD');
    % axis square;
    % 
    subplot(1,3,3);
    hold on;
    bar(1,dOrig,'FaceColor','r');
    bar(2,wSorig,'FaceColor','k');
    bar(3,rbThresholdOrig,'FaceColor',[0.5 0.5 0.5]);
    set(gca,'XTick',[1 2 3]);
    set(gca,'XTickLabel',{'d' 'D_{opt}' 'Threshold'});
    title('Weights');
    set(gca,'FontSize',20);
    ylim(max([dOrig wSorig rbThresholdOrig]).*[-1.2 1.2]);
    axis square;

    figure;
    set(gcf,'Position',[189 395 1280 420]);
    subplot(1,3,1);
    plot(c,meanChangeXvec,'ko','LineWidth',1);
    xlim([-3 3]);
    ylim([-3 3]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaD');
    ylabel('Measured \DeltaD');
    title(['Correlation = ' num2str(rhoColorSwitch,3)]);
    axis square;
    
    subplot(1,3,2);
    deltaApredictedNoColor = deltaS.*wS;
    hold on;
    for i = 1:length(deltaApredictedNoColor)
       if c(i)>0
           clrPlot = [0.6 0 0];
       elseif c(i)<0
           clrPlot = [0 0 1];
       else
           clrPlot = [0.6 0 1];
       end        
       plot(deltaApredictedNoColor(i),meanChangeXvec(i),'o','LineWidth',1,'Color',clrPlot,'MarkerFaceColor',clrPlot); 
    end
    xlim([-3 3]);
    ylim([-3 3]);
    set(gca,'FontSize',15);
    xlabel('Prediction \DeltaD');
    ylabel('Measured \DeltaD');
    axis square;
    
%     subplot(1,3,3);
%     hold on;
%     bar(1,wRorig,'FaceColor','r');
%     bar(2,wSorig,'FaceColor','k');
%     bar(3,biasOrig,'FaceColor','w');
%     bar(4,rbThresholdOrig,'FaceColor',[0.5 0.5 0.5]);
%     set(gca,'XTick',[1 2 3 4]);
%     set(gca,'XTickLabel',{'Red' 'D_{opt}' 'Bias' 'Threshold'});
%     title('Weights');
%     set(gca,'FontSize',20);
%     ylim(max([wRorig wSorig biasOrig rbThresholdOrig]).*[-1.2 1.2]);
%     axis square;    
end

d = dOrig;
wS = wSorig;
rbThreshold = rbThresholdOrig;

end