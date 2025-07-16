%%

bRandomize = false;
bPlot = true;
subjNumAll = [1 3 5 10 16 17 18 20];
% subjNumAllRand = subjNumAll;
peakLocModelPredictionAll = [];
nRepeat = 1;
corrShuffle = [];
subjNumInclude = [1:8];
fileStr = 'acuityModelingPrediction';

figure;
set(gcf,'Position',[176 273 1309 669]);
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
            subplot(2,4,j);
            hold on;
            plot(shiftValBestFit+defocusForStim+modelPrediction875nmPurpleAt2pt5-2.5,normcdf(dprimeMetric.*dprimeScaleBestFit/2),'-','Color',[0.56 0 1],'LineWidth',1);
            errorbar(unqFocDst.*scaleFac,normcdf(dprime/2),(normcdf(dprime/2)-normcdf(dprimeCI(1,:)/2)),(normcdf(dprimeCI(2,:)/2)-normcdf(dprime/2)),'o','Color',[0.56 0 1],'MarkerFaceColor','w','LineWidth',1.5,'MarkerSize',10);
            xlabel('Gabor - word Distance (D)');
            ylabel('Proportion Correct');
            set(gca,'FontSize',15);
            axis square;
            text(-1,0.3,['S' num2str(j)],'FontSize',18);
            % title(['Mean defocus at 875nm = ' num2str(modelPrediction875nmPurpleAt2pt5,3) 'D']);
            ylim([0.2 1]);
            xlim([-1.2 1.2]);
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
end
