%% RUN THIS SCRIPT AFTER ARCimgQualityAnalysisConesMeanZ

wvInd2examine = [21 26 31 36 41 46 51 56 61];
wave = 380:4:780;

figure;
set(gcf,'Position',[232 199 1199 722]);
for i = 1:length(wvInd2examine)
    subplot(3,3,i);
    imagesc(ifftshift(abs(ifft2(squeeze(oi.optics.OTF.OTF(:,:,wvInd2examine(i))))))); 
    axis square; 
    colormap gray;
    xlim([134 184]);
    ylim([81 130]);
    set(gca,'XTick',[]);
    set(gca,'YTick',[]);
    title(['\lambda = ' num2str(wave(wvInd2examine(i))) 'nm']);
end