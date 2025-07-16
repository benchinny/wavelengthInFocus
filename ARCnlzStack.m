function [x3stack,tInterp,AFCp,indCndCell] = ARCnlzStack(sn,bPLOT)

bEXCLUDE = true; % EXCLUDE BLINKS AND BAD TRIALS? 
bLeaveOutTransitions = true; % LEAVE OUT FIRST 50 FRAMES AND TRANSITION PERIOD OF ACCOMMODATION?

if sn==11 % 'VISIT' NUMBERS
   vs = [2 3 4 7];
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
   vs = 7:10;
   excludeTrials = [];     
elseif sn==17
   vs = 1:4;
   excludeTrials = [];    
elseif sn==18
   vs = 3:6;
   excludeTrials = [];       
elseif sn==19
   vs = 5:8;
   excludeTrials = [];    
elseif sn==21
   vs = 1:4;
   excludeTrials = [];  
elseif sn==23
   vs = 14:17;
   excludeTrials = [];   
elseif sn==24
   vs = 9:12;
   excludeTrials = [];      
elseif sn==26
   vs = 4:7;
   excludeTrials = [];    
elseif sn==25
   vs = 6:9;
   excludeTrials = [];       
elseif sn==27
   vs = 5:8;
   excludeTrials = [];
elseif sn==28
   vs = 13:16;
   excludeTrials = [];   
elseif sn==29
   vs = 7:10;
   excludeTrials = [];
elseif sn==30
   vs = 9:12;
   excludeTrials = [];         
elseif sn==31
   vs = 7:10;
   excludeTrials = [];      
elseif sn==32
   vs = 11:14;
   excludeTrials = []; 
elseif sn==33
   vs = [3 4 5 7];
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
elseif strcmp(getenv("USER"),'emily')
   dataDirectory = '/Users/emily/Library/CloudStorage/GoogleDrive-emilyacooper@gmail.com/Shared drives/ARChroma/Analysis/';
elseif strcmp(getenv("USER"),'ben')
   dataDirectory = '/home/ben/Documents/ARchroma/Analysis/';
elseif strcmp(getenv("USER"),'benjaminchin')
   dataDirectory = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Analysis/';   
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
       x3stack(indCnd(k0),1:length(x3interp)) = x3interp-x3interp(1);
       y3stack(indCnd(k0),1:length(y3tmp)) = y3tmp-y3tmp(1);
       % ACCOMMODATIVE DEMAND FROM EXPERIMENT 
       sinValuesTmp = AFCp.sinValues(indCnd(k0),:);
       tSin = 0:(1/(length(sinValuesTmp)-1)):1;
       tSinInterp = 0:(1/(length(x2{k0})-1)):1;
       accContinuous = interp1(tSin,sinValuesTmp,tSinInterp); 
       % THIS IS AN OBNOXIOUS WAY OF COMPUTING THE AVERAGE CHANGE WITHIN A
       % TRIAL, BUT IT WORKS AND MAY BE A BIT MORE ROBUST
       if bLeaveOutTransitions
          diffVec = imresize([0 -4/length(accContinuous) 0 4/length(accContinuous)],size(accContinuous),'nearest');
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

figSize = 3;

% SCRAMBLING CONDITIONS FOR ORGANIZING PLOTS LATER
indB2mixed = uniqueRGBvalues(:,1)<0.0001 & uniqueRGBvalues(:,2)<0.0001 ...
           & uniqueRGBvalues(:,3)>0.0001 & uniqueRGBvalues(:,4)>0.0001 ...
           & uniqueRGBvalues(:,5)<0.0001 & uniqueRGBvalues(:,6)>0.0001;

indR2mixed = uniqueRGBvalues(:,1)>0.0001 & uniqueRGBvalues(:,2)<0.0001 ...
           & uniqueRGBvalues(:,3)<0.0001 & uniqueRGBvalues(:,4)>0.0001 ...
           & uniqueRGBvalues(:,5)<0.0001 & uniqueRGBvalues(:,6)>0.0001;

indSame = (uniqueRGBvalues(:,1)>0.0001 & uniqueRGBvalues(:,2)<0.0001 ...
        & uniqueRGBvalues(:,3)<0.0001 & uniqueRGBvalues(:,4)>0.0001 ...
        & uniqueRGBvalues(:,5)<0.0001 & uniqueRGBvalues(:,6)<0.0001) ...
        | (uniqueRGBvalues(:,1)<0.0001 & uniqueRGBvalues(:,2)<0.0001 ...
        & uniqueRGBvalues(:,3)>0.0001 & uniqueRGBvalues(:,4)<0.0001 ...
        & uniqueRGBvalues(:,5)<0.0001 & uniqueRGBvalues(:,6)>0.0001);

indR2BorB2R = (uniqueRGBvalues(:,1)>0.0001 & uniqueRGBvalues(:,2)<0.0001 ...
             & uniqueRGBvalues(:,3)<0.0001 & uniqueRGBvalues(:,4)<0.0001 ...
             & uniqueRGBvalues(:,5)<0.0001 & uniqueRGBvalues(:,6)>0.0001) ...
            | (uniqueRGBvalues(:,1)<0.0001 & uniqueRGBvalues(:,2)<0.0001 ...
             & uniqueRGBvalues(:,3)>0.0001 & uniqueRGBvalues(:,4)>0.0001 ...
             & uniqueRGBvalues(:,5)<0.0001 & uniqueRGBvalues(:,6)<0.0001);

uniqueRGBvalues = [uniqueRGBvalues(indSame,:); uniqueRGBvalues(indR2BorB2R,:); uniqueRGBvalues(indB2mixed,:); uniqueRGBvalues(indR2mixed,:)];

if bPLOT
    for i = 1:size(uniqueRGBvalues,1)
        indRGB = ismember(uniqueConditions(:,1:6),uniqueRGBvalues(i,:),'rows');
        stepSizes = uniqueConditions(indRGB,7);
        rowNum = mod(i,figSize);
        if rowNum==0
           rowNum = figSize;
        end
        if rowNum==1
            figNum = floor(i/figSize)+1;
            figure(figNum);
            set(gcf,'Position',[207 189 1240 765]);
        end
        for j = 1:length(stepSizes)
            colNum = length(stepSizes)+1;
            indUnq = ismember(uniqueConditions(:,1:6),uniqueRGBvalues(i,:),'rows') ...
                      & abs(uniqueConditions(:,7)-stepSizes(j))<0.001;     
            indCndMarkers = indCndCell{indUnq};
            subplot(figSize,colNum,j+(rowNum-1)*colNum);
            set(gca,'FontSize',15);
            frmDuration = 0.033;
            hold on;
            plot([0:frmDuration:(size(x3stack,2)-1)*frmDuration],x3stack(indCndMarkers,:),'-','Color',[0 0.45 0.74]);
            % plot([0:frmDuration:(size(x3stack,2)-1)*frmDuration],y3stack(indCndMarkers,:),'-','Color',[0.85 0.33 0.1]);
            plot([0:frmDuration:(size(x3stack,2)-1)*frmDuration],mean(x3stack(indCndMarkers,:),1),'-','Color',[0 0.45 0.74],'LineWidth',2);
            % plot([0:frmDuration:(size(x3stack,2)-1)*frmDuration],mean(y3stack(indCndMarkers,:),1),'-','Color',[0.85 0.33 0.1],'LineWidth',2);
            % xlim([0 quantile(v1length,0.05)*frmDuration]);
            xlim([0 6]);
            ylim([-3 3]);
            xlabel('Time (s)'); ylabel('Relative Power (Diopters)'); 
            title(['Step = ' num2str(stepSizes(j)*optDistScale) ...
                  ', RGB = [' num2str(uniqueRGBvalues(i,1)) ' ' num2str(uniqueRGBvalues(i,2)) ' ' num2str(uniqueRGBvalues(i,3)) '] to ['...
                  num2str(uniqueRGBvalues(i,4)) ' ' num2str(uniqueRGBvalues(i,5)) ' ' num2str(uniqueRGBvalues(i,6)) ']']); 
        end
        subplot(figSize,length(stepSizes)+1,colNum+colNum*(rowNum-1));
        hold on;
        changeLabels = {};
        for j = 1:length(stepSizes)
            indUnq = ismember(uniqueConditions(:,1:6),uniqueRGBvalues(i,:),'rows') ...
                      & abs(uniqueConditions(:,7)-stepSizes(j))<0.001;   
            bar(j,mean(meanChangeX(indUnq,:)),'FaceColor','w');
            bar(j+length(stepSizes),mean(meanChangeY(indUnq,:)),'FaceColor','w');
            errorbar(j,mean(meanChangeX(indUnq,:)),std(meanChangeX(indUnq,:)),'k');
            errorbar(j+length(stepSizes),mean(meanChangeY(indUnq,:)),std(meanChangeY(indUnq,:)),'k');
            stepString = '-+';
            changeLabels{j} = ['H' stepString(int16(1+(1+sign(stepSizes(j)))/2))];
            changeLabels{j+length(stepSizes)} = ['V' stepString(int16(1+(1+sign(stepSizes(j)))/2))];
        end
        ylabel('Mean Accommodative Response (D)');
        set(gca,'FontSize',15);
        set(gca,'XTick',[1:length(stepSizes)*2]);
        set(gca,'XTickLabel',changeLabels);        
        ylim([-3 3]);
    end
end
