function AFCnlzSin         
sn=input('Subject Number? '); sn=sn+1000;
vs=input('Visit Number? ')
ey=input('Eye? 1forRight 2forBinocular ');
filePath = 'G:\My Drive\exp_bvams\code_repo\';

% LEAVE THE FOLLOWING TWO LINES OF CODE UNCOMMENTED UNTIL I FIGURE OUT 
% WHAT c0 IS FOR
% load('JnJrand211130', 'TCAm')
% c0=TCAm(sn-1000, vs-1);
if ey(1)==2
    load('G:\My Drive\exp_bvams\code_repo\ARC\AFCflsB.mat'); c1=2;
elseif ey(1)==1
    load('G:\My Drive\exp_bvams\code_repo\ARC\AFCflsR.mat'); c1=1;
end
fnm=AFCfls{sn-1000, vs}; load(fnm)
opticalDistanceScale = 0.82;

%  [jsonFile,jsonPath] = uigetfile('*.json','' , 'G:\.shortcut-targets-by-id\17-MjlIMJ6eySxBl-1ikM1jyh1UQiX26y\code_repo\JnJ')
[jsonFile,jsonPath] = uigetfile('*.json','' , 'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real') %stim conputer
%         [jsonFile,jsonPath] = uigetfile('*.json','' , 'G:\My Drive\m-wamvideocapture\videos\processed\alwaysrecord_notcp')

% dt=jsondecode(jsonFile);
dt=jsondecode(fileread([jsonPath jsonFile]));

%horizontal
x=-1.*[(cell2mat(struct2cell(dt.ext_right_mu))-cell2mat(struct2cell(dt.ext_left_mu)))];
y=-1.*[(cell2mat(struct2cell(dt.ext_bottom_mu))-cell2mat(struct2cell(dt.ext_top_mu)))];
           
t3=cell2mat(struct2cell(dt.time_totalsecs))
v0=[find(diff(t3)>1); length(t3)];
i0=1;
for k0=1:length(v0);
   v1=i0:v0(k0);
   x0{k0}=x(v1);
   y0{k0}=y(v1);
   i0=v0(k0)+1;
end

%i0=find(AFCv==3);
i0=find(AFCp.v00(:,1)==0);
i0 = 1;
           
x1=[]; y1=[];
for k0=1:length(i0)
 x1=[x1 x0{i0(k0)}'];
 y1=[y1 y0{i0(k0)}'];  
end
x1 = x1(1:10);
y1 = y1(1:10);
           
%x3=[]; y3=[]; v1=[-0.5:0.25:4]; v2=v1(AFCv); v3=[];
clear x2 y2
x3=[]; y3=[]; 
% v1=AFCp.v0; 
v2=AFCp.v00; 
meanv2 = AFCp.meanv00;
v3=[]; 
% x4=ones(1,length(v2)); x5=zeros(length(v1), 4); y5=x5; 

for k0=1:size(v2,1)
   x2{k0}=x0{k0}-mean(x1);
   y2{k0}=y0{k0}-mean(y1);
   x3=[x3 x2{k0}'];
   y3=[y3 y2{k0}'];
   sinValuesTmp = AFCp.sinValues(k0,:);
   tSin = 0:(1/(length(sinValuesTmp)-1)):1;
   tSinInterp = 0:(1/(length(x2{k0})-1)):1;
   v3 = [v3 opticalDistanceScale.*interp1(tSin,sinValuesTmp,tSinInterp)];
%   v3=[v3 v2(k0).*ones(1,length(x2{k0}))];
%    i1=find(v2(k0)==v1)
%    x5(i1, x4(i1))=mean(x2{k0}./-3.59);
%    y5(i1, x4(i1))=mean(y2{k0}./-3.33);
%    x6(k0)=mean(x2{k0}./-3.59);
%    y6(k0)=mean(y2{k0}./-3.33);
%    x4(i1)=x4(i1)+1;
end

x3=x3./-2.87; y3=y3./-2.58;
plot([1:length(x3)], [v3; x3; y3])
xlabel('Frame'); ylabel('Power (Diopters)'); title('Autorefractor measurement'); legend('Demand', 'Horizontal', 'Vertical')
           
% %figure; plot(v1, mean(x5,2))
% figure; errorbar(v1, mean(x5,2), std(x5,0, 2)./sqrt(4))
% hold on;
% errorbar(v1, mean(y5,2), std(y5,0, 2)./sqrt(4))
% xlabel('Demand'); ylabel('Response measured'); legend('Horizontal', 'Vertical'); title(fnm) 
