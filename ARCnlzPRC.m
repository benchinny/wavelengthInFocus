function ARCnlzPRC         
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
         
scaleFactor = 0.8;
 %        [jsonFile,jsonPath] = uigetfile('*.json','' , 'G:\.shortcut-targets-by-id\17-MjlIMJ6eySxBl-1ikM1jyh1UQiX26y\code_repo\JnJ')
         [jsonFile,jsonPath] = uigetfile('*.json','' , 'G:\My Drive\ARchromaVideoCapture\videos\processed\centralperipheral_real') %stim conputer
         %         [jsonFile,jsonPath] = uigetfile('*.json','' , 'G:\My Drive\m-wamvideocapture\videos\processed\alwaysrecord_notcp')

%         dt=jsondecode(jsonFile);
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
           i0=find(AFCp.meanv00==2.5);

           
           x1=[]; y1=[];
           for k0=1:length(i0)
             x1=[x1 x0{i0(k0)}'];
             y1=[y1 y0{i0(k0)}'];  
           end
           
           %x3=[]; y3=[]; v1=[-0.5:0.25:4]; v2=v1(AFCv); v3=[];
           clear x2 y2
           x3=[]; y3=[]; v1=AFCp.v0; v2=AFCp.meanv00; v3=[]; x4=ones(1,length(v2)); x5=zeros(length(v1), 4); y5=x5; 

           for k0=1:length(v2)
              x2{k0}=x0{k0}-mean(x1);
              y2{k0}=y0{k0}-mean(y1);
           x3=[x3 x2{k0}'];
           y3=[y3 y2{k0}'];
           v3=[v3 scaleFactor.*v2(k0).*ones(1,length(x2{k0}))];
           i1=find(v2(k0)==v1)
           x5(i1, x4(i1))=mean(x2{k0}./-3.59);
           y5(i1, x4(i1))=mean(y2{k0}./-3.33);
           x6(k0)=mean(x2{k0}./-3.59);
           y6(k0)=mean(y2{k0}./-3.33);
           x4(i1)=x4(i1)+1;
           
           end
           x3=x3./-3.59; y3=y3./-3.33;
           plot([1:length(x3)], [v3; x3; y3])
           xlabel('Frame'); ylabel('Power (Diopters)'); title('Autorefractor measurement'); legend('Demand', 'Horizontal', 'Vertical')
           
           
           %figure; plot(v1, mean(x5,2))
           figure; errorbar(v1, mean(x5,2), std(x5,0, 2)./sqrt(4))
           hold on;
           errorbar(v1, mean(y5,2), std(y5,0, 2)./sqrt(4))
           xlabel('Demand'); ylabel('Response measured'); legend('Horizontal', 'Vertical'); title(fnm) 

           v2' 
           x7=(x6+y6)./2
           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECTION COMMENTED OUT BY BEN FOR NOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% m=[repmat([sn-1000 vs c0 c1], [76 1]) [1:76]' v2' x7' x6' y6'];
% a=repmat(uint32(str2num(fnm(end-9:end))), [76 1]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUFF BELOW WAS ALREADY COMMENTED OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%m=[repmat([sn vs c0 c1], [8 1]) [1:8]' y0' v' vL' vH' repmat(double(a), [8 1])]
%m=[repmat([sn vs c0 c1], [76 1]) [1:76]' v2' x7' x6' y6'], [repmat(uint32(str2num(fnm(end-9:end))), [76 1])]]

% xlswrite('tstX.xlsx', m)
% m0=xlsread('Sample_VA.xlsx');
% xlswrite('tstX.xlsx', m, 'A10:I18')
%load('JnJ\ACMDflnm', 'ACMDflnm'); m0=xlsread(ACMDflnm); ACMDflnm=['JnJ\ACMDxls' fnm(end-9:end) '.xlsx'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMMENTED OUT BY BEN FOR NOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [afcFile,afcPath] = uigetfile('*.xlsx','' , 'G:\.shortcut-targets-by-id\17-MjlIMJ6eySxBl-1ikM1jyh1UQiX26y\code_repo\JnJ')
% m0=xlsread([afcPath afcFile]); ACMDflnm=['JnJ\ACMDxls' fnm(end-9:end) '.xlsx'];
% 
% xlswrite(ACMDflnm, m0);
% xlswrite(ACMDflnm, m, ['A' n2s(size(m0,1)+2) ':I' n2s(size(m0,1)+77)]);
% xlswrite(ACMDflnm, a, ['J' n2s(size(m0,1)+2) ':J' n2s(size(m0,1)+77)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STUFF BELOW WAS ALREADY COMMENTED OUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save('JnJ\ACMDflnm', 'ACMDflnm');

% 
% 
% m=[repmat([sn-1000 vs c0 c1], [8 1]) [1:8]' y0' v' vL' vH']
% a=repmat(uint32(str2num(fnm(end-9:end))), [8 1])
% 
% 
% [etmFile,etmPath] = uigetfile('*.xlsx','' , 'G:\.shortcut-targets-by-id\17-MjlIMJ6eySxBl-1ikM1jyh1UQiX26y\code_repo\JnJ')
% m0=xlsread([etmPath etmFile]); ETMflnm=['JnJ\ETMxls' fnm(end-9:end) '.xlsx'];
% 
% xlswrite(ETMflnm, m0);
% xlswrite(ETMflnm, m, ['A' n2s(size(m0,1)+2) ':I' n2s(size(m0,1)+9)])
% xlswrite(ETMflnm, a, ['J' n2s(size(m0,1)+2) ':J' n2s(size(m0,1)+9)])
