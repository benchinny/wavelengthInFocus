%%

rgb = [];
meanFocstmOptDst = [];
focStmOptDstIncr = [];
rspAcu = [];
stimOrientation = [];

% filePath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/ARChroma/Analysis/';
filePath = 'G:\My Drive\exp_bvams\code_repo\ARC\';
subj = 'BC';

if strcmp(subj,'IM')
    filenames = {[filePath 'S1005V16_AFC_RightACL0_2311141650.mat'] ...
                 [filePath 'S1005V15_AFC_RightACL0_2311141641.mat'] ...
                 [filePath 'S1005V14_AFC_RightACL0_2311141628.mat'] ...
                 [filePath 'S1005V13_AFC_RightACL0_2311141616.mat'] ...
                 [filePath 'S1005V12_AFC_RightACL0_2311141604.mat'] ...
                 [filePath 'S1005V18_AFC_RightACL0_2311151638.mat'] ...
                 [filePath 'S1005V19_AFC_RightACL0_2311151652.mat'] ...
                 [filePath 'S1005V20_AFC_RightACL0_2311151702.mat'] ...
                 [filePath 'S1005V21_AFC_RightACL0_2311151710.mat'] ...                 
                 [filePath 'S1005V22_AFC_RightACL0_2311151720.mat']};
end
if strcmp(subj,'BC')
%     filenames = {[filePath 'S1003V16_AFC_RightACL0_2311171554.mat'] ...
%                  [filePath 'S1003V16_AFC_RightACL0_2311171611.mat'] ...
%                  [filePath 'S1003V16_AFC_RightACL0_2311171624.mat'] ...
%                  [filePath 'S1003V16_AFC_RightACL0_2311171636.mat'] ...
%                  [filePath 'S1003V16_AFC_RightACL0_2311171730.mat'] ...
%                  [filePath 'S1003V16_AFC_RightACL0_2311171742.mat'] ...  
%                  [filePath 'S1003V19_AFC_RightACL0_2311251543.mat'] ...  
%                  [filePath 'S1003V19_AFC_RightACL0_2311251604.mat'] ...  
%                  [filePath 'S1003V19_AFC_RightACL0_2311251615.mat'] ...  
%                  [filePath 'S1003V19_AFC_RightACL0_2311251625.mat'] ...                   
%                  };
    filenames = {[filePath 'S1006V7_AFC_RightACL0_2312201030.mat'] ...
                 [filePath 'S1006V7_AFC_RightACL0_2312201044.mat'] ...   
                 [filePath 'S1006V7_AFC_RightACL0_2312201101.mat'] ...   
                 [filePath 'S1006V7_AFC_RightACL0_2312201129.mat'] ...   
                 [filePath 'S1006V7_AFC_RightACL0_2312201115.mat'] ...
                 [filePath 'S1006V7_AFC_RightACL0_2312201144.mat'] ...
                 };
    % filenames = {[filePath 'S1006V7_AFC_RightACL0_2312131147.mat'] ...
    %              [filePath 'S1006V7_AFC_RightACL0_2312131154.mat'] ...                  
    %              };
end

if strcmp(subj,'EC')
    filenames = {[filePath 'S1006V9_AFC_RightACL0_2312131507.mat'] ...
                 [filePath 'S1006V9_AFC_RightACL0_2312150944.mat'] ...
                 [filePath 'S1006V9_AFC_RightACL0_2312150952.mat'] ...
                 [filePath 'S1006V9_AFC_RightACL0_2312150958.mat'] ...
                 [filePath 'S1006V9_AFC_RightACL0_2312151005.mat'] ...
                 [filePath 'S1006V9_AFC_RightACL0_2312151011.mat'] ...  
                 [filePath 'S1006V9_AFC_RightACL0_2312151016.mat'] ...  
                 [filePath 'S1006V9_AFC_RightACL0_2312151021.mat'] ...   
                 [filePath 'S1006V9_AFC_RightACL0_2312181546.mat'] ...  
                 [filePath 'S1006V9_AFC_RightACL0_2312181553.mat'] ...  
                 [filePath 'S1006V9_AFC_RightACL0_2312181558.mat'] ...                    
                 };
end

for i = 1:length(filenames)
    load(filenames{i});
    if (i==3 || i==10) && strcmp(subj,'IM')
        rgb = [rgb; AFCp.rgb(1:131,:)];
        meanFocstmOptDst = [meanFocstmOptDst; AFCp.meanFocstmOptDst(1:131)];
        focStmOptDstIncr = [focStmOptDstIncr; AFCp.focStmOptDstIncr(1:131)];
        rspAcu = [rspAcu; AFCp.rspAcu(1:131)'];
        stimOrientation = [stimOrientation; AFCp.stimOrientation(1:131)];
    else
        rgb = [rgb; AFCp.rgb];
        meanFocstmOptDst = [meanFocstmOptDst; AFCp.meanFocstmOptDst];
        focStmOptDstIncr = [focStmOptDstIncr; AFCp.focStmOptDstIncr];
        rspAcu = [rspAcu; AFCp.rspAcu'];
        stimOrientation = [stimOrientation; AFCp.stimOrientation];
    end
end

meanFocstmOptDstUnq = unique(AFCp.meanFocstmOptDst);
indRed = rgb(:,1)>0 & rgb(:,3)==0;
indBlue = rgb(:,1)==0 & rgb(:,3)>0;
indMix = rgb(:,1)>0 & rgb(:,3)>0;
indRed1 = rgb(:,1)>0 & rgb(:,3)==0 & meanFocstmOptDst==meanFocstmOptDstUnq(1);
indBlue1 = rgb(:,1)==0 & rgb(:,3)>0 & meanFocstmOptDst==meanFocstmOptDstUnq(1);
indMix1 = rgb(:,1)>0 & rgb(:,3)>0 & meanFocstmOptDst==meanFocstmOptDstUnq(1);
if strcmp(subj,'BC') || strcmp(subj,'IM')
   indRed2 = rgb(:,1)>0 & rgb(:,3)==0 & meanFocstmOptDst==meanFocstmOptDstUnq(2);
   indBlue2 = rgb(:,1)==0 & rgb(:,3)>0 & meanFocstmOptDst==meanFocstmOptDstUnq(2);
   indMix2 = rgb(:,1)>0 & rgb(:,3)>0 & meanFocstmOptDst==meanFocstmOptDstUnq(2);
end
scaleFactor = 0.8;
focStmOptDstIncr = focStmOptDstIncr.*scaleFactor;

figure;
set(gcf,'Position',[248 499 1264 420]);
ind = indRed;
[~,~,~,~,PCdta,~,~] = psyfitgengauss(zeros(size(focStmOptDstIncr(ind))),focStmOptDstIncr(ind),rspAcu(ind)==stimOrientation(ind),[],[],[],1,4,0);
subplot(1,3,1);
plot(unique(focStmOptDstIncr(ind)),PCdta,'r-','LineWidth',1);
formatFigure('Relative distance (diopters)','Proportion Correct');
axis square;
ylim([0 1]);
ind = indBlue;
[~,~,~,~,PCdta,~,~] = psyfitgengauss(zeros(size(focStmOptDstIncr(ind))),focStmOptDstIncr(ind),rspAcu(ind)==stimOrientation(ind),[],[],[],1,4,0);
subplot(1,3,2);
plot(unique(focStmOptDstIncr(ind)),PCdta,'b-','LineWidth',1);
formatFigure('Relative distance (diopters)','Proportion Correct');
axis square;
ylim([0 1]);
ind = indMix;
[~,~,~,~,PCdta,~,~] = psyfitgengauss(zeros(size(focStmOptDstIncr(ind))),focStmOptDstIncr(ind),rspAcu(ind)==stimOrientation(ind),[],[],[],1,4,0);
subplot(1,3,3);
plot(unique(focStmOptDstIncr(ind)),PCdta,'m-','LineWidth',1);
formatFigure('Relative distance (diopters)','Proportion Correct');
axis square;
ylim([0 1]);

figure;
set(gcf,'Position',[248 499 1264 420]);
ind = indRed1;
[~,~,~,~,PCdta,~,~] = psyfitgengauss(zeros(size(focStmOptDstIncr(ind))),focStmOptDstIncr(ind),rspAcu(ind)==stimOrientation(ind),[],[],[],1,4,0);
subplot(1,3,1);
plot(unique(focStmOptDstIncr(ind)),PCdta,'r-','LineWidth',1);
formatFigure('Relative distance (diopters)','Proportion Correct');
axis square;
ylim([0 1]);
ind = indBlue1;
[~,~,~,~,PCdta,~,~] = psyfitgengauss(zeros(size(focStmOptDstIncr(ind))),focStmOptDstIncr(ind),rspAcu(ind)==stimOrientation(ind),[],[],[],1,4,0);
subplot(1,3,2);
plot(unique(focStmOptDstIncr(ind)),PCdta,'b-','LineWidth',1);
formatFigure('Relative distance (diopters)','Proportion Correct',['Accommodative distance = ' num2str(scaleFactor.*meanFocstmOptDstUnq(1)) 'D']);
axis square;
ylim([0 1]);
ind = indMix1;
[~,~,~,~,PCdta,~,~] = psyfitgengauss(zeros(size(focStmOptDstIncr(ind))),focStmOptDstIncr(ind),rspAcu(ind)==stimOrientation(ind),[],[],[],1,4,0);
subplot(1,3,3);
plot(unique(focStmOptDstIncr(ind)),PCdta,'m-','LineWidth',1);
formatFigure('Relative distance (diopters)','Proportion Correct');
axis square;
ylim([0 1]);

if strcmp(subj,'BC') || strcmp(subj,'IM')
    figure;
    set(gcf,'Position',[248 499 1264 420]);
    ind = indRed2;
    [~,~,~,~,PCdta,~,~] = psyfitgengauss(zeros(size(focStmOptDstIncr(ind))),focStmOptDstIncr(ind),rspAcu(ind)==stimOrientation(ind),[],[],[],1,4,0);
    subplot(1,3,1);
    plot(unique(focStmOptDstIncr(ind)),PCdta,'r-','LineWidth',1);
    formatFigure('Relative distance (diopters)','Proportion Correct');
    axis square;
    ylim([0 1]);
    ind = indBlue2;
    [~,~,~,~,PCdta,~,~] = psyfitgengauss(zeros(size(focStmOptDstIncr(ind))),focStmOptDstIncr(ind),rspAcu(ind)==stimOrientation(ind),[],[],[],1,4,0);
    subplot(1,3,2);
    plot(unique(focStmOptDstIncr(ind)),PCdta,'b-','LineWidth',1);
    formatFigure('Relative distance (diopters)','Proportion Correct',['Accommodative distance = ' num2str(scaleFactor.*meanFocstmOptDstUnq(2)) 'D']);
    axis square;
    ylim([0 1]);
    ind = indMix2;
    [~,~,~,~,PCdta,~,~] = psyfitgengauss(zeros(size(focStmOptDstIncr(ind))),focStmOptDstIncr(ind),rspAcu(ind)==stimOrientation(ind),[],[],[],1,4,0);
    subplot(1,3,3);
    plot(unique(focStmOptDstIncr(ind)),PCdta,'m-','LineWidth',1);
    formatFigure('Relative distance (diopters)','Proportion Correct');
    axis square;
    ylim([0 1]);
end
