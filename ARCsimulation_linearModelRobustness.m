%%

sn = 18; % CURRENTLY HAVE SUBJECTS 11 THROUGH 26
bEXCLUDE = true;
gammaFactorR = 2.4;
gammaFactorB = 2.4;
scaleFactor = 0.8;

if sn==11 % 'VISIT' NUMBERS
   vs = [2 3 4 5 6 7];
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
   vs = 7:12;
   excludeTrials = [];     
elseif sn==17
   vs = 1:6;
   excludeTrials = [];  
elseif sn==18
   vs = 3:10;
   excludeTrials = [];        
elseif sn==19
   vs = 5:10;
   excludeTrials = [];     
elseif sn==21
   vs = 1:6;
   excludeTrials = [];
elseif sn==23
   vs = 14:21;
   excludeTrials = [];   
elseif sn==24
   vs = 9:16;
   excludeTrials = [];     
elseif sn==25
   vs = 6:13;
   excludeTrials = [];        
elseif sn==26
   vs = 4:11;
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

deltaR = AFCp.rgb200(:,1).^gammaFactorR - AFCp.rgb100(:,1).^gammaFactorR;
deltaB = AFCp.rgb200(:,3).^gammaFactorB - AFCp.rgb100(:,3).^gammaFactorB    ;
deltaS = AFCp.v00*scaleFactor;
delta1 = ones(size(deltaR));

trials2use = 88:128;

deltaR = deltaR(trials2use);
deltaB = deltaB(trials2use);
deltaS = deltaS(trials2use);
delta1 = delta1(trials2use);

%%

nRepeats = 10000;

for i = 1:nRepeats
   testNoise = 0.6;
   weightsRBS1true = [0.6 -0.4 0.9];
   accRspSim = [deltaR deltaB deltaS]*weightsRBS1true' + normrnd(0,testNoise,size(deltaS));
   recoveredWeights(i,:) = [deltaR deltaB deltaS]\(accRspSim);
end

figure;
set(gcf,'Position',[407 246 920 645]);
subplot(2,2,1);
hist(recoveredWeights(:,1),101); hold on;
plot(weightsRBS1true(1).*[1 1],ylim,'k--','LineWidth',1);
xlabel('w_R');
subplot(2,2,2);
hist(recoveredWeights(:,2),101); hold on;
plot(weightsRBS1true(2).*[1 1],ylim,'k--','LineWidth',1);
xlabel('w_B');
subplot(2,2,3);
hist(recoveredWeights(:,3),101); hold on;
plot(weightsRBS1true(3).*[1 1],ylim,'k--','LineWidth',1);
xlabel('w_S');
