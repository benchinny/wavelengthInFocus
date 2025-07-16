function [defocus,strehl,timestamps] = ARCindividualAnalysisFIAT(subjNum)

if subjNum==1
    subjName = 'BenChin-OS';
    blockNums = [2 3 4 5 6];
    trialNums = [[1:20]' [1:20]' [1:20]' [1:20]' [1:20]'];
    % blockNums = [2 3];
    % trialNums = [[1:20]' [1:20]'];    
end

blockNum = randsample(blockNums,1);
trialNum = randsample(trialNums(:,blockNums==blockNum),1);

AFCp = ARCloadFileBVAMS(subjNum,blockNum);
[defocus, strehl, timestamps] = ARCloadFileFIAT(subjName,blockNum,trialNum,1);

timestamps = double(seconds(timestamps));
defocus = -defocus;
strehl = -strehl;

figure; 
hold on;
plot(timestamps-timestamps(1),imresize(AFCp.sinValues(trialNum,:)./1.25,[1 length(defocus)],'nearest'),'-k');
plot(timestamps-timestamps(1),defocus,'-k');
plot(timestamps-timestamps(1),strehl,'--k');
ylim([1 4]);
formatFigure('Time (s)','Defocus (D)',['Block ' num2str(blockNum) ' Trial ' num2str(trialNum)]);
plot(1.5,1.25,'s','MarkerFaceColor',AFCp.rgb100(trialNum,:),'MarkerSize',25,'Color',AFCp.rgb100(trialNum,:));
plot(4.5,1.25,'s','MarkerFaceColor',AFCp.rgb200(trialNum,:),'MarkerSize',25,'Color',AFCp.rgb200(trialNum,:));

end