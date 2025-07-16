function [defocus,timestamps] = ARCacuAnalysisFIAT(subjNum)

if subjNum==1
    subjName = 'BenChin-OS';
    blockNums = 7;
end

defocus = {};
timestamps = {};

for i = 1:length(blockNums)
    AFCp = ARCloadFileBVAMS(subjNum,blockNums(i));
    for j = 1:size(AFCp.t3,1)
        [defocusTmp, ~, timestampsTmp] = ARCloadFileFIAT(subjName,blockNums(i),j,0);
        defocus{end+1} = defocusTmp;
        timestamps{end+1} = double(seconds(timestampsTmp));
    end
end

focStmOptDstIncrUnq = unique(AFCp.focStmOptDstIncr);

figure;
set(gcf,'Position',[414 303 959 628]);
for i = 1:length(focStmOptDstIncrUnq)
    ind = find(abs(AFCp.focStmOptDstIncr-focStmOptDstIncrUnq(i))<0.001);
    subplot(3,3,i);
    hold on;
    for j = 1:length(ind)
        timestampsTrial = timestamps{ind(j)};
        defocusTrial = -defocus{ind(j)};
        indBad = defocusTrial==0;
        plot(timestampsTrial(~indBad)-timestampsTrial(1),defocusTrial(~indBad),'k');
    end
    ylim([2 4]);
    formatFigure('Time(s)','Defocus',['Increment ' num2str(focStmOptDstIncrUnq(i)./1.25) 'D']);
end

end