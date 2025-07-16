%%

bRandomize = false;
bPlot = true;
subjNumAll = [1 3 5 10 16 17 18 20];
% subjNumAllRand = subjNumAll;
peakLocModelPredictionAll = [];
nRepeat = 1;
corrShuffle = [];
subjNumInclude = [1:8];
fileStr = 'acuityModelingPredictionDiffLim';

for k = 1:nRepeat
    for j = 1:length(subjNumAll)
    
        subjNum = subjNumAll(j);
        
        load(['/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/acuityModeling/' fileStr 'S' num2str(subjNum) '.mat']);
        
        if bRandomize
            subjNumAllRand = subjNumAll(randperm(length(subjNumAll)));
            load(['/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/acuityModeling/' fileStr 'S' ...
                  num2str(subjNumAllRand(j)) '.mat'],'dprimeMetric','defocusForStim');
        end
    
        scaleFac = 0.8;
        % dprime2regress = interp1(defocusForStim+modelPrediction875nmPurpleAt2pt5,dprimeMetric,2.5+unqFocDst.*scaleFac);
        % 
        % % dprimeScale = max(dprime(:)./max(dprimeMetric));
        % dprimeScale = dprime2regress\dprime';
        
        shiftVals = -0.25:0.05:0.25;
        
        for i = 1:length(shiftVals)
            dprime2regressTmp = interp1(defocusForStim+modelPrediction875nmPurpleAt2pt5+shiftVals(i),dprimeMetric,2.5+unqFocDst.*scaleFac);
            dprimeScaleTmp(i) = dprime2regressTmp\dprime';
            errorDP(i) = sqrt(mean((dprimeScaleTmp(i).*dprime2regressTmp-dprime').^2));
        end
        [~,indMinShift] = min(errorDP);
        shiftValBestFit = shiftVals(indMinShift);
        dprimeScaleBestFit = dprimeScaleTmp(indMinShift);
        
        if bPlot
            figure;
            set(gcf,'Position',[342 460 1052 440]);
            subplot(1,2,1);
            hold on;
            plot(shiftValBestFit+defocusForStim+modelPrediction875nmPurpleAt2pt5,dprimeMetric.*dprimeScaleBestFit,'-','Color',[0.56 0 1],'LineWidth',1);
            errorbar(2.5+unqFocDst.*scaleFac,dprime,(dprime-dprimeCI(1,:)),(dprimeCI(2,:)-dprime),'o','Color',[0.56 0 1],'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
            xlabel('Distance');
            ylabel('D-prime metric');
            set(gca,'FontSize',15);
            axis square;
            title(['Mean defocus at 875nm = ' num2str(modelPrediction875nmPurpleAt2pt5,3) 'D']);
            % subplot(1,2,2);
            % plot(wvInFocusForStim,dprimeMetric,'k-','LineWidth',1);
            % xlabel('Wavelength in focus (nm)');
            % ylabel('D-prime metric');
            % set(gca,'FontSize',15);
            % axis square;
            
            figure;
            set(gcf,'Position',[342 460 1052 440]);
            subplot(1,2,1);
            hold on;
            plot(shiftValBestFit+defocusForStim+modelPrediction875nmPurpleAt2pt5-2.5,normcdf(dprimeMetric.*dprimeScaleBestFit/2),'-','Color',[0.56 0 1],'LineWidth',1);
            xlabel('Distance');
            ylabel('Proportion Correct');
            set(gca,'FontSize',15);
            axis square;
            title(['Mean defocus at 875nm = ' num2str(modelPrediction875nmPurpleAt2pt5,3) 'D']);
            ylim([0.4 1]);
            % subplot(1,2,2);
            % plot(wvInFocusForStim,dprimeMetric,'k-','LineWidth',1);
            % xlabel('Wavelength in focus (nm)');
            % ylabel('D-prime metric');
            % set(gca,'FontSize',15);
            % axis square;
        end
        
        stimDistanceSmp = 1.2:0.01:3.8;
        dprimeMetricSmooth = interp1(defocusForStim+modelPrediction875nmPurpleAt2pt5+shiftVals(indMinShift),dprimeMetric,stimDistanceSmp,'spline');
        [~,indPeak] = max(dprimeMetricSmooth);
        peakLocModelPrediction = stimDistanceSmp(indPeak);
        peakLocModelPredictionAll(j) = peakLocModelPrediction;
    end
    
    %%
    
    % peakLocModelPredictionAll = [1.84 2.07 2.39 2.19 2.23 1.92 1.97 2.14];
    peakLocActualAll = [1.8635 1.9435 2.6335 2.1235 2.6435 1.9035 1.8935 2.3235];
    peakLocLBall = [1.6835 1.8135 1.5935 1.9435 2.2535 1.7935 1.5035 2.2635];
    peakLocUBall = [3.0935 3.1835 3.5635 3.1035 2.7735 2.0535 2.4535 2.4235];
    
    if mod(k,20)==0
        figure;
        hold on;
        plot(peakLocModelPredictionAll(subjNumInclude),peakLocActualAll(subjNumInclude),'ko','MarkerFaceColor',[0.56 0 1],'MarkerSize',15);
        % errorbar(peakLocModelPredictionAll(subjNumInclude), ...
        %          peakLocActualAll(subjNumInclude), ...
        %          peakLocActualAll(subjNumInclude)-peakLocLBall(subjNumInclude), ...
        %          peakLocUBall(subjNumInclude)-peakLocActualAll(subjNumInclude), ...
        %          'ko','MarkerFaceColor',[0.56 0 1],'MarkerSize',15);
        axis square;
        set(gca,'FontSize',15);
        xlim([1.5 2.8]);
        ylim([1.5 2.8]);
        plot([1.5 2.8],[1.5 2.8],'k--','LineWidth',1);
        xlabel(['Predicted Peak Location (D)']);
        ylabel(['Actual Peak Location (D)']);
        title(['Correlation = ' num2str(corr(peakLocModelPredictionAll(subjNumInclude)',peakLocActualAll(subjNumInclude)'),3)]);
    end
    corrShuffle(k) = corr(peakLocModelPredictionAll',peakLocActualAll');
    display(['Iteration ' num2str(k)]);
end

%%

load('/Users/benjaminchin/Documents/ARchromaScraps/corrShuffle.mat','corrShuffle');

edgeShuffle = linspace(-0.2,1,25);
hShuffle = histogram(corrShuffle,edgeShuffle','FaceColor','w','LineWidth',1);
set(gca,'FontSize',15);
xlabel('Correlation');
hold on;
plot([0.886 0.886],ylim,'k--','LineWidth',1);

%%

peakLocModelPredictionAll = [1.8400 2.0700 2.3900 2.1900 2.2300 1.9200 1.9700 2.1400];
peakLocActualAll = [1.8635 1.9435 2.6335 2.1235 2.6435 1.9035 1.8935 2.3235];

dataTableAcuity = array2table([peakLocModelPredictionAll' peakLocActualAll'], ...
                             'VariableNames',{'Predictions' 'Actual'});

