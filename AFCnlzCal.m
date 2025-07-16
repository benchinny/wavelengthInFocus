%%

sn=1; sn=sn+1000;
vs=7;
ey=1;
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

fileNamesAll = {'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real\2023-03-24 16.03.json' ...
                'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real\2023-03-24 16.15.json' ...
                'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real\2023-03-24 16.17.json' ...
                'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real\2023-03-24 16.20.json' ...
                'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real\2023-03-24 16.24.json' ...
                'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real\2023-03-24 16.27.json' ...
                'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real\2023-03-24 16.29.json'};

x = [];
y = [];
t3 = [];
for i = 1:length(fileNamesAll)
    % dt=jsondecode(jsonFile);
    dt=jsondecode(fileread(fileNamesAll{i}));
    %horizontal
    xTmp=1.*[(cell2mat(struct2cell(dt.ext_right_mu))-cell2mat(struct2cell(dt.ext_left_mu)))];
    yTmp=1.*[(cell2mat(struct2cell(dt.ext_bottom_mu))-cell2mat(struct2cell(dt.ext_top_mu)))]; 
    t3tmp=cell2mat(struct2cell(dt.time_totalsecs));
    x = [x; xTmp];
    y = [y; yTmp];
    t3 = [t3; t3tmp];
end

% AFCp.v00 = [100/3.34 - 100./(0.1.*(43.53-[0:0.5:7])) ...
%             100/3.34 - 100./(0.1.*(35.95-[0:0.5:7 7.5 8.0 8.5 8.5 8.5]))];
AFCp.v00 = [(1.*(43.53-[0:0.5:7])) (1.*(35.95-[0:0.5:7 7.5 8.0 8.5 8.5 8.5]))];
% AFCp.v00 = [100/3.34 - 100./(0.1.*(43.53-[0:0.5:9])) ...
%             100/3.34 - 100./(0.1.*(35.95-[0:0.5:7.5]))];
% AFCp.v00 = [1:35];
anchorValue = AFCp.v00(1);
AFCp.v00 = AFCp.v00-anchorValue;

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
           
x1=[]; y1=[];
for k0=1:length(i0)
 x1=[x1 x0{i0(k0)}'];
 y1=[y1 y0{i0(k0)}'];  
end
           
%x3=[]; y3=[]; v1=[-0.5:0.25:4]; v2=v1(AFCv); v3=[];
clear x2 y2
x3=[]; y3=[]; 
% v1=AFCp.v0; 
v2=AFCp.v00; 
meanv2 = AFCp.meanv00;
v3=[]; 
% x4=ones(1,length(v2)); x5=zeros(length(v1), 4); y5=x5; 

for k0=1:length(v2)
%   x2{k0}=x0{k0}-mean(x1);
%   y2{k0}=y0{k0}-mean(y1);
   x2{k0}=x0{k0};
   y2{k0}=y0{k0};   
   x3=[x3 x2{k0}'];
   y3=[y3 y2{k0}'];
%   sinValuesTmp = AFCp.sinValues(k0,:);
%   tSin = 0:(1/(length(sinValuesTmp)-1)):1;
%   tSinInterp = 0:(1/(length(x2{k0})-1)):1;
%   v3 = [v3 interp1(tSin,sinValuesTmp,tSinInterp)];
   v3=[v3 v2(k0).*ones(1,length(x2{k0}))];
%    i1=find(v2(k0)==v1)
%    x5(i1, x4(i1))=mean(x2{k0}./-3.59);
%    y5(i1, x4(i1))=mean(y2{k0}./-3.33);
%    x6(k0)=mean(x2{k0}./-3.59);
%    y6(k0)=mean(y2{k0}./-3.33);
%    x4(i1)=x4(i1)+1;
end

x3=x3./1; y3=y3./1;
plot([1:length(x3)], [v3; x3; y3])
xlabel('Frame'); ylabel(''); title('Autorefractor measurement'); legend('Demand (D)', 'Right - left bars (pix)', 'Bottom - top bars (pix)')

figure;
hold on;
for i = 1:length(x2)
   plot(AFCp.v00(i)+anchorValue,mean(x2{i}),'ro');
   plot(AFCp.v00(i)+anchorValue,mean(y2{i}),'bo');
end
set(gca,'FontSize',15);
xlabel('Distance from lens to retina (mm)'); ylabel('Distance between bars (pixels)');
legend('Left - Right bars','Top - bottom bars');
% %figure; plot(v1, mean(x5,2))
% figure; errorbar(v1, mean(x5,2), std(x5,0, 2)./sqrt(4))
% hold on;
% errorbar(v1, mean(y5,2), std(y5,0, 2)./sqrt(4))
% xlabel('Demand'); ylabel('Response measured'); legend('Horizontal', 'Vertical'); title(fnm) 
