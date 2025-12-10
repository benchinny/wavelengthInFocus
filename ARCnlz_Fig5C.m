%% FIGURE 5C IN MANUSCRIPT

dprimeRatioAll = [2.41 1.000 1.000 1.41 1.000 1.98 1.76 1.023];

figure;
boxplot(dprimeRatioAll);
set(gca,'FontSize',15);
hold on;
for i = 1:length(dprimeRatioAll)
    if abs(dprimeRatioAll(i)-1)<0.01
        plot(1+i*0.1,dprimeRatioAll(i),'ko','MarkerSize',10,'MarkerFaceColor',[0.56 0 1]);
    else
        plot(1,dprimeRatioAll(i),'ko','MarkerSize',10,'MarkerFaceColor',[0.56 0 1]);
    end
end
ylabel('d''_{max}/d''_{0}');
