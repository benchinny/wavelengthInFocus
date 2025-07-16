%% SCRIPT FOR ANALYZING LCA DATA OF ALL PARTICIPANTS WHO PASSED SCREENING

clear;

%%

% subjNumAll = [1 3 5 9 10 16 17 18 20];
subjNumAll = [1 3 5 10 16 17 18 20];
defocusLCAmeasuredAll = [];
q1bestAll = [];
q2bestAll = [];
q3bestAll = [];
D0all = [];
defocusLCAmeasuredBootsAll = [];

for i = 1:length(subjNumAll)
    [defocusLCAmeasured, q1best, q2best, q3best,defocusLCAmeasuredBoots,D0] = ARCacuAnalysisLCA(subjNumAll(i),0,100);
    defocusLCAmeasuredAll(i,:) = defocusLCAmeasured;
    q1bestAll(i) = q1best;
    q2bestAll(i) = q2best;
    q3bestAll(i) = q3best;
    defocusLCAmeasuredBootsAll(:,:,i) = defocusLCAmeasuredBoots;
    D0all(i) = D0;
    display(['Finished subject ' num2str(i)]);
end

%%

wave2fit = 460:5:875;
standardLCAfit = humanWaveDefocus(wave2fit);
val532 = humanWaveDefocus(532);
standardLCAfit = standardLCAfit-val532;
displacementLambda = -4:4;
wavePlot = 380:5:875;
defocusAt875centered550 = [];

figure;
hold on;
% plot(wave2fit,standardLCAfit,'k-','LineWidth',1);
for i = 1:size(defocusLCAmeasuredAll,1)
    set(gca,'ColorOrderIndex',i);
    defocusAt875centered550(i) = humanWaveDefocusARC(533,wavePlot(end),subjNumAll(i)) + D0all(i)+defocusLCAmeasuredAll(i,2);
    plot(wavePlot,(-2-defocusAt875centered550(i))+humanWaveDefocusARC(533,wavePlot,subjNumAll(i))+D0all(i)+defocusLCAmeasuredAll(i,2),'-','LineWidth',1.5);
    defocus875to620(i,:) = humanWaveDefocusARC(533,wavePlot([48 100]),subjNumAll(i))+D0all(i)+defocusLCAmeasuredAll(i,2);
end
for i = 1:size(defocusLCAmeasuredAll,1)
    defocusLCAmeasuredBootsTmp = squeeze(defocusLCAmeasuredBootsAll(:,:,i))-1*defocusLCAmeasuredAll(i,2);
    CIlca = quantile(defocusLCAmeasuredBootsTmp',[0.025 0.975]);
    defocusLCAmeasuredAll(i,:) = defocusLCAmeasuredAll(i,:);
    set(gca,'ColorOrderIndex',i);
    errorbar([616 533 468]-displacementLambda(i), ...
             -(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2))+(-2-defocusAt875centered550(i)), ...
             -(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2)-CIlca(2,:)), ...
             -(CIlca(1,:)-(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2))), ...
             '.','MarkerSize',10, ...
        'LineWidth',1.5,'MarkerFaceColor','w');    
end
for i = 1:size(defocusLCAmeasuredAll,1)
    defocusLCAmeasuredBootsTmp = squeeze(defocusLCAmeasuredBootsAll(:,:,i));
    CIlca = quantile(defocusLCAmeasuredBootsTmp',[0.025 0.975]);   
    set(gca,'ColorOrderIndex',i);
    plot([616 533 468]-displacementLambda(i),-(defocusLCAmeasuredAll(i,:)-defocusLCAmeasuredAll(i,2))+(-2-defocusAt875centered550(i)),'o','MarkerSize',10, ...
        'LineWidth',1.5,'MarkerFaceColor','w');
end
% axis square;
xlim([440 875]);
% ylim([-1.5 2]);
formatFigure('Wavelength (\lambda)','Defocus (D)');
legend('','','','','','','','','','','','','','','','','S1','S2','S3','S4','S5','S6','S7','S8','Location','SouthEast');

figure;
hold on;
plot(ones([1 8]),defocus875to620(:,2)-defocus875to620(:,1),'ko','LineWidth',1);
plot(1,humanWaveDefocus(875)-humanWaveDefocus(620),'k*','LineWidth',1);
plot(1,mean(defocus875to620(:,2)-defocus875to620(:,1)),'kp','LineWidth',1);
ylim([0 0.8]);
ylabel('D_{875}-D_{620}');
set(gca,'FontSize',15);
set(gca,'XTick',[]);
