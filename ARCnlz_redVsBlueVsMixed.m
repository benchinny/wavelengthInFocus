function [x3stackRed,x3stackBlue,x3stackMixed,x3stackMixedMoreRed,x3stackMixedMoreBlue,meanRedPerTrial,meanBluePerTrial,meanMixedPerTrial,meanMixedMoreRedPerTrial,meanMixedMoreBluePerTrial,tInterp,AFCp,indCndCell] = ARCnlz_redVsBlueVsMixed(sn,bPLOT,lumCutoff)

bEXCLUDE = true; % EXCLUDE BLINKS AND BAD TRIALS? 
bLeaveOutTransitions = true; % LEAVE OUT FIRST 50 FRAMES AND TRANSITION PERIOD OF ACCOMMODATION?

if sn==11 % 'VISIT' NUMBERS
   vs = [2 3 4 5 6 7];
   excludeTrials = [];
elseif sn==12
   vs = [3:7];
   excludeTrials = [87 7 3 4 5 90 91 10 6 86 88];
elseif sn==13
   vs = 1:4;
   excludeTrials = [59 39 38 58 77 56 60 74];
elseif sn==14
   vs = 1:4;
   excludeTrials = [];    
elseif sn==15
   vs = [7:10];
   excludeTrials = [49];  
elseif sn==16
   vs = 7:12;
   excludeTrials = [];     
elseif sn==17
   vs = 1:4;
   excludeTrials = [];    
elseif sn==18
   vs = 3:10;
   excludeTrials = [];       
elseif sn==19
   vs = 5:8;
   excludeTrials = [];    
elseif sn==21
   vs = 1:4;
   excludeTrials = [];  
elseif sn==23
   vs = 14:21;
   excludeTrials = [];   
elseif sn==24
   vs = 9:16;
   excludeTrials = [];      
elseif sn==26
   vs = 4:11;
   excludeTrials = [];    
elseif sn==25
   vs = 6:9;
   excludeTrials = [];       
elseif sn==27
   vs = 5:12;
   excludeTrials = [];
elseif sn==28
   vs = 13:20;
   excludeTrials = [];   
elseif sn==29
   vs = 7:14;
   excludeTrials = [];   
elseif sn==30
   % vs = 4:11;
   vs = 9:16;
   excludeTrials = [];        
elseif sn==32
   vs = 11:18;
   excludeTrials = []; 
elseif sn==33
   vs = [3 4 5 7:11];
   excludeTrials = []; 
else
   error('ARCnlz: unhandled subject number!');
end
ey = 1; % 1 for Right eye, 2 for Binocular ');
filePath = 'G:\My Drive\exp_bvams\code_repo\';


if strcmp(getenv('username'),'bankslab')
   dataDirectory = [filePath 'ARC\'];
elseif strcmp(getenv("USER"),'benchin')
   dataDirectory = '/Users/benchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Analysis/'; 
elseif strcmp(getenv("USER"),'benjaminchin')
   dataDirectory = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Analysis/';
elseif strcmp(getenv("USER"),'emily')
   dataDirectory = '/Users/emily/Library/CloudStorage/GoogleDrive-emilyacooper@gmail.com/Shared drives/ARChroma/Analysis/';
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
       AFCpAll = AFCp;
    else
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
v1length = [];
for k0=1:length(v0);
   v1=i0:v0(k0);
   x0{k0}=x(v1);
   y0{k0}=y(v1);
   t0{k0}=t3(v1);
   i0=v0(k0)+1;
   v1length(end+1) = length(v1);
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
x1 = x1(11:21);
y1 = y1(11:21);

% MAIN ANALYSIS
uniqueConditions = unique([AFCp.rgb100 AFCp.rgb200 AFCp.v00],'rows');
uniqueRGBvalues = unique([AFCp.rgb100 AFCp.rgb200],'rows');
xScale = -2.87; % SCALE FACTOR CONVERTING PIXELS TO DIOPTERS
yScale = -2.58; % SCALE FACTOR CONVERTING PIXELS TO DIOPTERS
optDistScale = 0.8;
rgbValuesAll = {};
x3stack = zeros([length(AFCp.AFCv) 400]);
y3stack = zeros([length(AFCp.AFCv) 400]);
tInterp = 0:0.033:6.6;

for i = 1:size(uniqueConditions,1)
    % INDEX OF COLOR CONDITIONS
    indCnd = find(ismember([AFCp.rgb100 AFCp.rgb200 AFCp.v00],uniqueConditions(i,:),'rows'));
    if bEXCLUDE
       indCnd = indCnd(~ismember(indCnd,excludeTrials));
    end
    rgbValues = [];
    for k0=1:size(indCnd,1)
       % ANALYZING SUBJECT'S ACCOMMODATION
       x2{k0}=x0{indCnd(k0)}-mean(x1); % MEAN CENTERING
       y2{k0}=y0{indCnd(k0)}-mean(y1); % MEAN CENTERING
       t2{k0}=t0{indCnd(k0)};
       x3tmp = (x2{k0}')./xScale;
       y3tmp = (y2{k0}')./yScale;
       % REMOVING OUTLIERS
       x3diff = [0 diff(x3tmp)];
       y3diff = [0 diff(y3tmp)];
       x3outliers = abs(x3diff)>0.7 | abs(x3tmp)>5;
       y3outliers = abs(y3diff)>0.7 | abs(y3tmp)>5;
       meanx3 = mean(x3tmp(~x3outliers));
       meany3 = mean(y3tmp(~y3outliers));
       if bEXCLUDE
           x3tmp(x3outliers) = meanx3;
           y3tmp(y3outliers) = meany3;
       end       
       x3tmp = [x3tmp 0];
       y3tmp = [y3tmp 0];
       t3tmp = t2{k0}';
       x3interp = interp1([t3tmp-t3tmp(1) max(tInterp)],[x3tmp],tInterp);
       x3stack(indCnd(k0),1:length(x3interp)) = x3interp;
       y3stack(indCnd(k0),1:length(y3tmp)) = y3tmp;
       % ACCOMMODATIVE DEMAND FROM EXPERIMENT 
       sinValuesTmp = AFCp.sinValues(indCnd(k0),:);
       tSin = 0:(1/(length(sinValuesTmp)-1)):1;
       tSinInterp = 0:(1/(length(x2{k0})-1)):1;
       accContinuous = interp1(tSin,sinValuesTmp,tSinInterp); 
       % THIS IS AN OBNOXIOUS WAY OF COMPUTING THE AVERAGE CHANGE WITHIN A
       % TRIAL, BUT IT WORKS AND MAY BE A BIT MORE ROBUST
       if bLeaveOutTransitions
          diffVec = imresize([0 4/length(accContinuous) 0 0],size(accContinuous),'nearest');
       else
           diffVec = imresize([-2/length(accContinuous) 2/length(accContinuous)],size(accContinuous),'nearest');
           if abs(corr(accContinuous',diffVec'))<0.95
               error('ARCnlz: you may want to check whether the step change occurs halfway through the trial, or not!');
           end
       end
       meanChangeX(i,k0) = sum(diffVec.*x2{k0}')./xScale;
       meanChangeY(i,k0) = sum(diffVec.*y2{k0}')./yScale;
       % VECTOR OF RGB VALUES FOR PLOTTING
       rgbValues = [rgbValues imresize([AFCp.rgb100(indCnd(k0),:)' AFCp.rgb200(indCnd(k0),:)'],[3 length(tSinInterp)],'nearest')];
    end

    % STORE ACCOMMODATION AND COLOR VALUES FOR PLOTTING
    rgbValuesAll{i} = rgbValues;
    indCndCell{i} = indCnd;
end

whichInterval = 1;
scaleEquateRB = 4;
gammaFactorR = 2.4;
gammaFactorB = 2.2;
maxLumCdm2 = 0.87;

if whichInterval==1
    rgb00 = AFCp.rgb100;
elseif whichInterval==2
    rgb00 = AFCp.rgb200;
end

unqRedValues = unique(rgb00(:,1));
indRedOnly = abs(rgb00(:,1)-max(unqRedValues))<0.001 & rgb00(:,2)==0 & rgb00(:,3)==0;
indBlueOnly = rgb00(:,1)==0 & rgb00(:,2)==0 & rgb00(:,3)>0;
indMixed = rgb00(:,1)>0 & rgb00(:,2)==0 & rgb00(:,3)>0 & ...
           maxLumCdm2.*(scaleEquateRB.*rgb00(:,1).^gammaFactorR + rgb00(:,3).^gammaFactorB)>lumCutoff(1) & ...
           maxLumCdm2.*(scaleEquateRB.*rgb00(:,1).^gammaFactorR + rgb00(:,3).^gammaFactorB)<lumCutoff(2) & ...
           rgb00(:,1)<0.62 & rgb00(:,3)>0.57;
indMixedMoreRed = rgb00(:,1)>0 & rgb00(:,2)==0 & rgb00(:,3)>0 & ...
                  scaleEquateRB.*rgb00(:,1).^gammaFactorR > rgb00(:,3).^gammaFactorB;
indMixedMoreBlue = rgb00(:,1)>0 & rgb00(:,2)==0 & rgb00(:,3)>0 & ...
                  scaleEquateRB.*rgb00(:,1).^gammaFactorR < rgb00(:,3).^gammaFactorB;

lengthHalf = 87;
if whichInterval==1
    x3stackRed = x3stack(indRedOnly,1:lengthHalf);
    x3stackBlue = x3stack(indBlueOnly,1:lengthHalf);
    x3stackMixed = x3stack(indMixed,1:lengthHalf);
    x3stackMixedMoreRed = x3stack(indMixedMoreRed,1:lengthHalf);
    x3stackMixedMoreBlue = x3stack(indMixedMoreBlue,1:lengthHalf);
elseif whichInterval==2
    x3stackRed = x3stack(indRedOnly,(lengthHalf+6):end);
    x3stackBlue = x3stack(indBlueOnly,(lengthHalf+6):end);
    x3stackMixed = x3stack(indMixed,(lengthHalf+6):end); 
    x3stackMixedMoreRed = x3stack(indMixedMoreRed,(lengthHalf+6):end);
    x3stackMixedMoreBlue = x3stack(indMixedMoreBlue,(lengthHalf+6):end);
end

meanMatrixRed = imresize([0 2],size(x3stackRed),'nearest');
meanMatrixBlue = imresize([0 2],size(x3stackBlue),'nearest');
meanMatrixMixed = imresize([0 2],size(x3stackMixed),'nearest');
meanMatrixMixedMoreRed = imresize([0 2],size(x3stackMixedMoreRed),'nearest');
meanMatrixMixedMoreBlue = imresize([0 2],size(x3stackMixedMoreBlue),'nearest');
meanRedPerTrial = mean(meanMatrixRed.*x3stackRed,2);
meanBluePerTrial = mean(meanMatrixBlue.*x3stackBlue,2);
meanMixedPerTrial = mean(meanMatrixMixed.*x3stackMixed,2);
meanMixedMoreRedPerTrial = mean(meanMatrixMixedMoreRed.*x3stackMixedMoreRed,2);
meanMixedMoreBluePerTrial = mean(meanMatrixMixedMoreBlue.*x3stackMixedMoreBlue,2);

if bPLOT
    figure;
    set(gca,'FontSize',15);
    frmDuration = 0.033;
    hold on;
    meanAnchor = mean(x3stack(indRedOnly,1));
%    plot([0:frmDuration:(lengthHalf-1)*frmDuration],x3stack(indRedOnly,1:lengthHalf),'-','Color',[1 0 0]);
    plot([0:frmDuration:(lengthHalf-1)*frmDuration],mean(x3stackRed)-meanAnchor,'-','Color',[1 0 0],'LineWidth',2);
    xlim([0 3]);
    ylim([-3 3]);
    xlabel('Time (s)'); ylabel('Relative Power (Diopters)');     
    axis square;
    plot([0:frmDuration:(lengthHalf-1)*frmDuration],mean(x3stackBlue)-meanAnchor,'-','Color',[0 0 1],'LineWidth',2);
%    plot([0:frmDuration:(lengthHalf-1)*frmDuration],x3stack(indMixed,1:lengthHalf),'-','Color',[1 0 1]);
    plot([0:frmDuration:(lengthHalf-1)*frmDuration],mean(x3stackMixed)-meanAnchor,'-','Color',[1 0 1],'LineWidth',2);
end

end