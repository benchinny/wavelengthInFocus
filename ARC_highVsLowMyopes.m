%%

low_myopes = [0.139 -0.288; ...
              0.218 -0.118; ...
              0.237 -0.433];

high_myopes = [0.713 -0.504; ...
               0.293 -0.544; ...
               0.419 -0.644; ...
               0.359 -0.384; ...
               0.201 -0.181; ...
               0.134 0.013; ...
              -0.143 -0.384];

meanHighMyopes = mean(high_myopes(:,1)+high_myopes(:,2));
stdErrorHighMyopes = std(high_myopes(:,1)+high_myopes(:,2))./sqrt(size(high_myopes,1));

meanLowMyopes = mean(low_myopes(:,1)+low_myopes(:,2));
stdErrorLowMyopes = std(low_myopes(:,1)+low_myopes(:,2))./sqrt(size(low_myopes,1));

figure;
hold on;
bar(1,meanHighMyopes,'w','LineWidth',1);
bar(2,meanLowMyopes,'w','LineWidth',1);
errorbar(1,meanHighMyopes,stdErrorHighMyopes,'k','LineWidth',1);
errorbar(2,meanLowMyopes,stdErrorLowMyopes,'k','LineWidth',1);
axis square;
set(gca,'FontSize',15);
set(gca,'XTick',[1 2])
set(gca,'XTickLabel',{'High myopes' 'Low myopes'});
ylabel('w_{R}-|w_{B}|');
