%%

subjNum = 10;
if subjNum==10
    subjName = 'S20-OD';
end
blockNumAll = 3:8;
trialNumAll = 1:36;
wLunq = -1:0.2:1;
wMunq = -1:0.2:1;
wSunq = [-0.6 0.2];
RMSE = [];
defocus875all = [];
defocus875predAll = [];
wvInFocusAll = [];
savePath = '/Users/benjaminchin/Library/CloudStorage/GoogleDrive-bechin@berkeley.edu/Shared drives/CIVO_BVAMS/data/coneWeightsError/';

for Sindex = 1:length(wSunq)
    for i = 1:length(wLunq)
        for j = 1:length(wMunq)
            wL = wLunq(i);
            wM = wMunq(j);
            wS = wSunq(Sindex);
            defocus875 = zeros([length(blockNumAll) length(trialNumAll)]);
            defocus875pred = zeros([length(blockNumAll) length(trialNumAll)]); 
            for k = 1:length(blockNumAll)
                parfor l = 1:length(trialNumAll)
                    wvInFocus(k,l) = ARCwvInFocusCones(subjNum,blockNumAll(k),trialNumAll(l),[wL wM wS]);
                    display(['wLind = ' num2str(wLunq(i)) ', wMind = ' num2str(wMunq(j)) ', wSind = ' num2str(wSunq(Sindex)) ', block number ' num2str(blockNumAll(k)) ', trial number ' num2str(trialNumAll(l))]);  
                    % LOAD ZERNIKE TABLE AND TIMESTAMPS
                    [ZernikeTable, ~, ~, TimeStamp] = ARCloadFileFIAT(subjName,blockNumAll(k),trialNumAll(l),0);
            
                    NumCoeffs = width(ZernikeTable)-8; % determine how many coefficients are in the cvs file. 
                    c=zeros(size(ZernikeTable,1),65); %this is the vector that contains the Zernike polynomial coefficients. We can work with up to 65. 
                    PARAMS = struct;
                    PARAMS.PupilSize=mean(table2array(ZernikeTable(:,5))); %default setting is the pupil size that the Zernike coeffs define, PARAMS(3)
                    PARAMS.PupilFitSize=mean(table2array(ZernikeTable(:,5))); 
                    PARAMS.PupilFieldSize=PARAMS.PupilSize*2; %automatically compute the field size
                    c(:,3:NumCoeffs)=table2array(ZernikeTable(:,11:width(ZernikeTable)));
                    indBad = c(:,4)==0;
                    meanC = mean(c(~indBad,:),1); % TAKE MEAN OF COEFFICIENTS    
                    defocusCorrectionFactor = (1e6/(4*sqrt(3)))*((PARAMS.PupilSize/2000)^2);
                    defocus875(k,l) = meanC(4)./defocusCorrectionFactor;
                    AFCp = ARCloadFileBVAMS(subjNum+10,blockNumAll(k));
                    defocus875pred(k,l) = AFCp.meanv00(trialNumAll(l))./1.2255 - humanWaveDefocusS10(wvInFocus(k,l),875);
                end
            end
            RMSE(i,j,Sindex) = sqrt(mean((defocus875pred(:)-defocus875(:)).^2));
            defocus875predAll(:,:,i,j,Sindex) = defocus875pred;
            defocus875all(:,:,i,j,Sindex) = defocus875;
            wvInFocusAll(:,:,i,j,Sindex) = wvInFocus;
        end
    end
   save([savePath 'wvInFocusModelResults' num2str(round(wSunq(Sindex)*10)) '.mat'],'defocus875all','defocus875predAll','wvInFocusAll','RMSE','wSunq'); 
end

% %%
% 
% figure; 
% hold on;
% histogram(wvInFocus(1,:),11,'FaceColor','r');
% histogram(wvInFocus(2,:),11,'FaceColor','g');
% histogram(wvInFocus(3,:),11,'FaceColor','b');
