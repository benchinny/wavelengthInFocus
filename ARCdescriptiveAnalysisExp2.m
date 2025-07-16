%% LOADING DATA FOR ORIGINAL ANALYSIS

subjNum = 1;

if subjNum==1
    subjName = 'BenChin-OS';
    blockNums = [2 3 4 5 6];
    trialNums = {[1:20]' [1:20]' [1:20]' [1:20]' [ 1:20]'};
    % blockNums = [2 3];
    % trialNums = [[1:20]' [1:20]']; 
    limVals = [-2 2];
elseif subjNum==2
    subjName = 'S2-OS';
    blockNums = [2 3 4 5 6];
    trialNums = {[1:20]' [1:20]' [1:20]' [1:20]' [1:20]'};
    % blockNums = [2 3];
    % trialNums = [[1:20]' [1:20]'];     
    limVals = [-2 2];
end

wvInFocus1all = [];
meanv00all = [];
v00all = [];
rgb1all = [];
rgb2all = [];
defocusBasic = [];
defocusBasic2 = [];
defocus875stack = [];
% STIMULUS DISTANCE
meanv00stack = [];
% STIMULUS COLOR
rgb1stack = [];
% CORRECTION FACTOR THAT IS ALSO IN AUSTIN'S CODE
defocusCorrectionFactor = 0.57735;

for l = 1:5 % LOOP OVER BLOCK
    trialsTmp = trialNums{l};
    for k = 1:length(trialsTmp) % LOOP OVER TRIAL
        % LOADING DATA
        blockNumInd = l;
        blockNumTmp = blockNums(blockNumInd);
        trialNumTmp = trialsTmp(k);
        
        AFCp = ARCloadFileBVAMS(subjNum,blockNumTmp); % LOAD BVAMS DATA
        % LOAD ZERNIKE TABLE AND TIMESTAMPS
        [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumTmp,trialNumTmp,0);
        % GET THE TIMESTAMP CORRESPONDING TO THE HALFWAY POINT
        t = seconds(TimeStamp)-min(seconds(TimeStamp));
        tHalfway = max(t)/2;
        tDiffFromHalfway = abs(t-tHalfway);
        [~,indMinT] = min(tDiffFromHalfway);
        FrameStart = (indMinT-29):indMinT; % analyze 30 frames
        FrameEnd = (length(t)-29):length(t);

        NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
        c=zeros(30,65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS.PupilSize=mean(table2array(ZernikeTable(FrameStart,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS.PupilFitSize=mean(table2array(ZernikeTable(FrameStart,5))); 
        PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
        c(:,3:NumCoeffs)=table2array(ZernikeTable(FrameStart,11:width(ZernikeTable)));
        indBad = c(:,4)==0;
        c(indBad,4) = mean(c(~indBad,4));
        meanC = mean(c,1); % TAKE MEAN OF COEFFICIENTS
        
        c2=zeros(30,65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
        PARAMS2.PupilSize=mean(table2array(ZernikeTable(FrameEnd,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
        PARAMS2.PupilFitSize=mean(table2array(ZernikeTable(FrameEnd,5))); 
        PARAMS2.PupilFieldSize=PARAMS2.PupilSize*2; %automatically compute the field size
        c2(:,3:NumCoeffs)=table2array(ZernikeTable(FrameEnd,11:width(ZernikeTable)));
        indBad = c2(:,4)==0;
        c2(indBad,4) = mean(c2(~indBad,4));        
        meanC2 = mean(c2,1); % TAKE MEAN OF COEFFICIENTS   

        % figure;
        % subplot(1,1,1);
        % plot(c2(:,4));
        % title(['Block ' num2str(blockNumTmp) ', Trial ' num2str(trialNumTmp)]);
        % ylim([0 2]);
        % pause;
        % close;

        % STORE COLORS FOR FIRST AND SECOND STIMULI
        rgb00 = [];
        rgb00(1,:) = AFCp.rgb100(trialNumTmp,:);
        rgb00(2,:) = AFCp.rgb200(trialNumTmp,:);
        defocusBasic(end+1,:) = meanC(4);
        defocusBasic2(end+1,:) = meanC2(4);
        rgb1all(end+1,:) = AFCp.rgb100(trialNumTmp,:);
        rgb2all(end+1,:) = AFCp.rgb200(trialNumTmp,:);
        meanv00all(end+1,:) = AFCp.meanv00(trialNumTmp);
        v00all(end+1,:) = AFCp.v00(trialNumTmp);

        defocus875stack(end+1,:) = meanC(4);
        meanv00stack(end+1,:) = AFCp.meanv00(trialNumTmp);
        rgb1stack(end+1,:) = AFCp.rgb100(trialNumTmp,:);
    end
end

% ORIGINAL ANALYSIS

scaleEquateRB = 4.1086;
scaleEquateRG = 9.6592;
gammaFactorR = 2.4;
gammaFactorG = 2.6;
gammaFactorB = 2.2;
maxLumCdm2 = 0.40;

deltaR = scaleEquateRB.*rgb2all(:,1).^gammaFactorR - scaleEquateRB.*rgb1all(:,1).^gammaFactorR;
deltaG = scaleEquateRG.*rgb2all(:,2).^gammaFactorG - scaleEquateRG.*rgb1all(:,2).^gammaFactorG;
deltaB = rgb2all(:,3).^gammaFactorB - rgb1all(:,3).^gammaFactorB;
% COMPONENTS OF LINEAR REGRESSION
deltaS = v00all*0.87;
delta1 = ones(size(deltaR));
deltaR = deltaR.*maxLumCdm2;
deltaB = deltaB.*maxLumCdm2;
deltaG = deltaG.*maxLumCdm2;
defocusChange = defocusBasic2./defocusCorrectionFactor - defocusBasic./defocusCorrectionFactor;
deltaR = deltaR(1:100);
deltaB = deltaB(1:100);
deltaS = deltaS(1:100);
defocusChange = defocusChange(1:100);
weightsRBS1 = [deltaR deltaB deltaS]\(defocusChange);

% COMPUTING CORRELATIONS BETWEEN DIFFERENT MODEL PREDICTIONS AND DATA
rhoFull = corr([deltaR deltaB deltaS]*weightsRBS1,defocusChange);
rhoNoColor = corr([deltaS]*weightsRBS1(3),defocusChange);
rhoColor = corr([deltaR deltaB]*weightsRBS1(1:2),defocusChange);

nParams = 4;
nParamsNoColor = 2;

% COMPUTING COMPONENETS FOR AIC
trialMeans = [deltaR deltaB deltaS]*weightsRBS1;
errorIndividual = defocusChange - trialMeans;
for i = 1:100
   [stdTmp(i),LLtmp(i)] = ARCfitStdGauss(errorIndividual);
end
[~,bestInd] = min(LLtmp);
estResidualStd = stdTmp(bestInd);
LL = sum(log(normpdf(defocusChange,trialMeans,estResidualStd)));
aic = 2*nParams-2*LL;

% COMPUTING COMPONENTS FOR AIC FOR NO-COLOR MODEL
trialMeansNoColor = [deltaS]*weightsRBS1(3);
errorIndividualNoColor = defocusChange - trialMeansNoColor;
for i = 1:100
   [stdTmp(i),LLtmp(i)] = ARCfitStdGauss(errorIndividualNoColor);
end
[~,bestInd] = min(LLtmp);
estResidualStdNoColor = stdTmp(bestInd);
LLnoColor = sum(log(normpdf(defocusChange,trialMeansNoColor,estResidualStdNoColor)));
aicNoColor = 2*nParamsNoColor-2*LLnoColor;

figure;
set(gcf,'Position',[189 395 1280 420]);
subplot(1,2,1);
%    plot([deltaR deltaB deltaS]*weightsRBS1_x,defocusChange,'ko','LineWidth',1);
predictionTmp = [deltaR deltaB deltaS]*weightsRBS1;
hold on;
for i = 1:length(predictionTmp)
    deltaRtmp = (1/maxLumCdm2)*deltaR(i);
    deltaBtmp = (1/maxLumCdm2)*deltaB(i);
    RBratio = 0.25.*(deltaRtmp-deltaBtmp+2);
    RBratio
    plot(predictionTmp(i),defocusChange(i),'o','LineWidth',1,'Color',[RBratio 0 1-RBratio],'MarkerFaceColor',[RBratio 0 1-RBratio]);
    % plot(predictionTmp(i),defocusChange(i),'o','LineWidth',1,'Color','k','MarkerFaceColor','w');
end    
plot(limVals,limVals,'k--');
xlim(limVals);
ylim(limVals);
set(gca,'FontSize',15);
xlabel('Prediction \DeltaA');
ylabel('Measured \DeltaA');
title(['Correlation = ' num2str(rhoFull,3)]);
axis square;

subplot(1,2,2);
hold on;
bar(1,weightsRBS1(1),'FaceColor','r');
bar(2,weightsRBS1(2),'FaceColor','b');
bar(3,weightsRBS1(3),'FaceColor','k');
set(gca,'XTick',[1 2 3]);
set(gca,'XTickLabel',{'Red' 'Blue' 'D_{opt}'});
title('Weights');
set(gca,'FontSize',20);
ylim(max(weightsRBS1).*[-1.2 1.2]);
axis square;

figure;
set(gcf,'Position',[189 395 1280 420]);
subplot(1,2,1);
% plot([deltaS]*weightsS1_x,defocusChange,'ko','LineWidth',1);
predictionTmp = [deltaS]*weightsRBS1(3);
hold on;
for i = 1:length(predictionTmp)
    deltaRtmp = (1/maxLumCdm2)*deltaR(i);
    deltaBtmp = (1/maxLumCdm2)*deltaB(i);
    RBratio = 0.25.*(deltaRtmp-deltaBtmp+2);
    RBratio
    plot(predictionTmp(i),defocusChange(i),'o','LineWidth',1,'Color',[RBratio 0 1-RBratio],'MarkerFaceColor',[RBratio 0 1-RBratio]);
    % plot(predictionTmp(i),defocusChange(i),'o','LineWidth',1,'Color','k','MarkerFaceColor','w');
end
plot(limVals,limVals,'k--');
xlim(limVals);
ylim(limVals);
set(gca,'FontSize',15);
xlabel('Prediction \DeltaA');
ylabel('Measured \DeltaA');
title(['Correlation = ' num2str(rhoNoColor,3)]);
axis square;

subplot(1,2,2);
hold on;
bar(1,weightsRBS1(3),'FaceColor','k');
set(gca,'XTick',[1 2]);
set(gca,'XTickLabel',{'D_{opt}'});
title('Weights');
set(gca,'FontSize',20);
ylim(max(weightsRBS1(3)).*[-1.2 1.2]);
axis square;

deltaSunq = unique(deltaS);

figure;
set(gcf,'Position',[214 355 1345 603]);
for i = 1:length(deltaSunq)
    subplot(2,4,i);
    hold on;
    ind = abs(deltaS-deltaSunq(i))<0.001;
    weightsRB1 = [deltaR(ind) deltaB(ind) ones([sum(ind) 1])]\(defocusChange(ind));
    weightsRB1all(:,i) = weightsRB1;
    predictionTmp = [deltaR(ind) deltaB(ind) ones([sum(ind) 1])]*weightsRB1;
    defocusChangeTmp = defocusChange(ind);
    rhoColor(i) = corr(predictionTmp,defocusChangeTmp);
    for j = 1:length(predictionTmp)
        deltaRtmp = (1/maxLumCdm2)*deltaR(j);
        deltaBtmp = (1/maxLumCdm2)*deltaB(j);
        RBratio = 0.25.*(deltaRtmp-deltaBtmp+2);
        plot(predictionTmp(j),defocusChangeTmp(j),'o','LineWidth',1,'Color',[RBratio 0 1-RBratio],'MarkerFaceColor',[RBratio 0 1-RBratio]);
        % plot(predictionTmp(j),defocusChange(j),'o','LineWidth',1,'Color','k','MarkerFaceColor','w');
    end    
    set(gca,'FontSize',12);
    xlabel('Prediction');
    ylabel('Measured defocus');
    title(['Correlation = ' num2str(rhoColor(i),3)]);
end

figure;
set(gcf,'Position',[214 355 1345 603]);
for i = 1:size(weightsRB1all,2)
    subplot(2,4,i);
    hold on;
    bar(1,weightsRB1all(1,i),'FaceColor','r');
    bar(2,weightsRB1all(2,i),'FaceColor','b');
    bar(3,weightsRB1all(3,i),'FaceColor','k');
    set(gca,'XTick',[1 2 3]);
    set(gca,'XTickLabel',{'Red' 'Blue' 'D_{opt}'});
    title('Weights');
    set(gca,'FontSize',20);
    ylim(max(abs(weightsRB1all(:,i))).*[-1.2 1.2]);
    axis square;
end

%% CORRELATION ANALYSES

rLum1 = 0.4.*scaleEquateRB.*rgb1stack(:,1).^2.4;
gLum1 = 0.4.*scaleEquateRG.*rgb1stack(:,2).^2.6;
bLum1 = 0.4.*rgb1stack(:,3).^2.2;
opticalDistCorrected = meanv00stack.*0.87-0.46;
opticalDistCorrectedUnq = unique(opticalDistCorrected);

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(rLum1(ind),defocus875stack(ind),'ro','MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(rLum1(ind),defocus875stack(ind)))]);      
    title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(gLum1(ind),defocus875stack(ind),'go','MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(gLum1(ind),defocus875stack(ind)))]);      
    title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(bLum1(ind),defocus875stack(ind),'bo','MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(bLum1(ind),defocus875stack(ind)))]);      
    title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(rLum1(ind)./bLum1(ind),defocus875stack(ind),'o','Color',[1 0 0.5],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(rLum1(ind)./bLum1(ind),defocus875stack(ind)))]);      
    if i==1
       title(['R/B, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(bLum1(ind)./rLum1(ind),defocus875stack(ind),'o','Color',[0.5 0 1],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(bLum1(ind)./rLum1(ind),defocus875stack(ind)))]);      
    if i==1
       title(['B/R, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(rLum1(ind)./gLum1(ind),defocus875stack(ind),'o','Color',[1 0.5 0],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(rLum1(ind)./gLum1(ind),defocus875stack(ind)))]);       
    if i==1
       title(['R/G, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(gLum1(ind)./rLum1(ind),defocus875stack(ind),'o','Color',[0.5 0.7 0],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(gLum1(ind)./rLum1(ind),defocus875stack(ind)))]);       
    if i==1
       title(['G/R, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(gLum1(ind)./bLum1(ind),defocus875stack(ind),'o','Color',[0 1 0.5],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(gLum1(ind)./bLum1(ind),defocus875stack(ind)))]);     
    if i==1
       title(['G/B, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(bLum1(ind)./gLum1(ind),defocus875stack(ind),'o','Color',[0 0.5 1],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(bLum1(ind)./gLum1(ind),defocus875stack(ind)))]);    
    if i==1
       title(['B/G, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(rLum1(ind)+bLum1(ind)+gLum1(ind),defocus875stack(ind),'o','Color',[0 0 0],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(rLum1(ind)+bLum1(ind)+gLum1(ind),defocus875stack(ind)))]);
    if i==1
       title(['R+G+B, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(rLum1(ind)-bLum1(ind),defocus875stack(ind),'o','Color',[1 0 0.5],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(rLum1(ind)-bLum1(ind),defocus875stack(ind)))]);      
    if i==1
       title(['R-B, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(rLum1(ind)-gLum1(ind),defocus875stack(ind),'o','Color',[1 0.5 0],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(rLum1(ind)-gLum1(ind),defocus875stack(ind)))]);      
    if i==1
       title(['R-G, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end

figure;
set(gcf,'Position',[267 401 1096 527]);
for i = 1:length(opticalDistCorrectedUnq)
    ind = abs(opticalDistCorrected-opticalDistCorrectedUnq(i))<0.001;
    subplot(2,3,i);
    plot(gLum1(ind)-bLum1(ind),defocus875stack(ind),'o','Color',[0 1 0.5],'MarkerSize',12,'MarkerFaceColor','w','LineWidth',1);
    axis square;
    set(gca,'FontSize',12);
    xlabel('Red luminance');
    ylabel('Defocus at 875nm');
    yLims = ylim;
    xLims = xlim; 
    text(xLims(1)+(xLims(2)-xLims(1))*0.5,yLims(1)+(yLims(2)-yLims(1))*0.75,['\rho = ' num2str(corr(gLum1(ind)-bLum1(ind),defocus875stack(ind)))]);      
    if i==1
       title(['G-B, Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    else
       title(['Optical distance = ' num2str(opticalDistCorrectedUnq(i),3)]);
    end
end
