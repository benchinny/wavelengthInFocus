%%

sn = [18 23 24 26 27 28 29 30 32 33];
xStackAll = [];
rgb100all = [];
rgb200all = [];
v00all = [];

for i = 1:length(sn)
    [x3stack,tInterp,AFCp] = ARCnlzStack(sn(i),0);
    xStackAll = [xStackAll; x3stack];
    rgb100all = [rgb100all; AFCp.rgb100];
    rgb200all = [rgb200all; AFCp.rgb200];
    v00all = [v00all; AFCp.v00];
end

uniqueConditions = unique([rgb100all rgb200all v00all],'rows');
uniqueRGBvalues = unique([rgb100all rgb200all],'rows');

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
optDistScale = 0.8;

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
        indCndMarkers = abs(rgb100all(:,1)-uniqueRGBvalues(i,1))<0.0001 & ...
                        abs(rgb100all(:,2)-uniqueRGBvalues(i,2))<0.0001 & ...
                        abs(rgb100all(:,3)-uniqueRGBvalues(i,3))<0.0001 & ...
                        abs(rgb200all(:,1)-uniqueRGBvalues(i,4))<0.0001 & ...
                        abs(rgb200all(:,2)-uniqueRGBvalues(i,5))<0.0001 & ...
                        abs(rgb200all(:,3)-uniqueRGBvalues(i,6))<0.0001 & ...
                        abs(v00all(:,1)-stepSizes(j))<0.0001;
        subplot(figSize,colNum,j+(rowNum-1)*colNum);
        set(gca,'FontSize',15);
        frmDuration = 0.033;
        hold on;
        % plot([0:frmDuration:(size(xStackAll,2)-1)*frmDuration],xStackAll(indCndMarkers,:),'-','Color',[0 0.45 0.74]);
        ind1stHalf = 1:91;
        ind2ndHalf = 91:182;
        tSamples = [0:frmDuration:(size(xStackAll,2)-1)*frmDuration];
        xStackCI1 = quantile(xStackAll(indCndMarkers,ind1stHalf),[0.16 0.84]);
        xStackCI2 = quantile(xStackAll(indCndMarkers,ind2ndHalf),[0.16 0.84]);
        fill([tSamples(ind1stHalf) fliplr(tSamples(ind1stHalf))],[xStackCI1(1,:) fliplr(xStackCI1(2,:))],uniqueRGBvalues(i,1:3),'EdgeColor','none');
        plot(tSamples(ind1stHalf),mean(xStackAll(indCndMarkers,ind1stHalf),1),'-','Color',uniqueRGBvalues(i,1:3),'LineWidth',2);
        fill([tSamples(ind2ndHalf) fliplr(tSamples(ind2ndHalf))],[xStackCI2(1,:) fliplr(xStackCI2(2,:))],uniqueRGBvalues(i,4:6),'EdgeColor','none');
        plot(tSamples(ind2ndHalf),mean(xStackAll(indCndMarkers,ind2ndHalf),1),'-','Color',uniqueRGBvalues(i,4:6),'LineWidth',2);
        xlim([0 6]);
        ylim([-0.75 2.8]);
        xlabel('Time (s)'); ylabel('Relative Power (Diopters)'); 
        title(['Step = ' num2str(stepSizes(j)*optDistScale) ...
              ', RGB = [' num2str(uniqueRGBvalues(i,1)) ' ' num2str(uniqueRGBvalues(i,2)) ' ' num2str(uniqueRGBvalues(i,3)) '] to ['...
              num2str(uniqueRGBvalues(i,4)) ' ' num2str(uniqueRGBvalues(i,5)) ' ' num2str(uniqueRGBvalues(i,6)) ']']); 
    end
%     subplot(figSize,length(stepSizes)+1,colNum+colNum*(rowNum-1));
%     hold on;
%     changeLabels = {};
%     for j = 1:length(stepSizes)
%         indUnq = ismember(uniqueConditions(:,1:6),uniqueRGBvalues(i,:),'rows') ...
%                   & abs(uniqueConditions(:,7)-stepSizes(j))<0.001;   
%         bar(j,mean(meanChangeX(indUnq,:)),'FaceColor','w');
%         bar(j+length(stepSizes),mean(meanChangeY(indUnq,:)),'FaceColor','w');
%         errorbar(j,mean(meanChangeX(indUnq,:)),std(meanChangeX(indUnq,:)),'k');
%         errorbar(j+length(stepSizes),mean(meanChangeY(indUnq,:)),std(meanChangeY(indUnq,:)),'k');
%         stepString = '-+';
%         changeLabels{j} = ['H' stepString(int16(1+(1+sign(stepSizes(j)))/2))];
%         changeLabels{j+length(stepSizes)} = ['V' stepString(int16(1+(1+sign(stepSizes(j)))/2))];
%     end
%     ylabel('Mean Accommodative Response (D)');
%     set(gca,'FontSize',15);
%     set(gca,'XTick',[1:length(stepSizes)*2]);
%     set(gca,'XTickLabel',changeLabels);        
%     ylim([-3 3]);
end
