%% SETUP PARAMETERS

subjNumAll = [1 3 5 10 16 17 18 20];
% subjNumAll = [3 10];
stimNum = 8;
LplusMweights = 0.0:0.1:1.0;
Sweights = fliplr(LplusMweights);
coneWeightsFolder = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsErrorSpatFilter/colorMechPredictions/';

%% MAKE PLOTS SHOWING IMAGE QUALITY AS A FUNCTION OF WAVELENGTH

for i = 1:length(subjNumAll)
    figure; 
    set(gcf,'Position',[278 182 1213 781]);
    for j = 1:length(LplusMweights)
        [wvInFocus, wave, peakCorr] = ARCwvInFocusConesMeanZspatFilter(subjNumAll(i),stimNum,[LplusMweights(j).*[1 1] -1*Sweights(j)]);
        display(['Subject ' num2str(subjNumAll(i)) ', L+M weight ' num2str(LplusMweights(j)) ', -S weight ' num2str(Sweights(j))]);
        subplot(3,4,j);
        plot(wave,peakCorr,'k-','LineWidth',1); hold on;
        plot([wvInFocus wvInFocus],[min(peakCorr)-0.1 max(peakCorr)+0.1],'k--');
        ylim([min(peakCorr)-0.1 max(peakCorr)+0.1]);
        axis square;
        set(gca,'FontSize',12);
        title(['L+M weight ' num2str(LplusMweights(j)) ', -S weight ' num2str(Sweights(j))]);
        if j==1
            ylabel('Peak correlation');
            title(['Subject ' num2str(subjNumAll(i))]);
        end
        if j==9
            xlabel('Wavelength in focus (\lambda)');
        end
    end
    saveas(gcf,[coneWeightsFolder 'LplusMminusSpred' num2str(subjNumAll(i))],'png');
end
