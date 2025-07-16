function [RMSE, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCtestWvInFocusMeanZspatFilterPlotHelper(subjNum,defocus875,rgbAll,optDistAll,w)

bFitMeans = true;

rgbUnq = unique(rgbAll,'rows');
wvInFocus = zeros([size(rgbUnq,1) 1]);
for l = 1:size(rgbUnq,1)
    wvInFocus(l,:) = ARCwvInFocusConesMeanZspatFilter(subjNum,l,w);
    if wvInFocus(l,:)<400 || wvInFocus(l,:)>720
        wvInFocus(l,:) = NaN;
    end
    display(['stim ' num2str(l)]);
end 

if ~bFitMeans
    wvInFocusAll = zeros([length(optDistAll) 1]);
    defocus875predTmp = zeros([length(optDistAll) 1]);
    for l = 1:size(rgbUnq,1)
        indStiml = abs(rgbAll(:,1)-rgbUnq(l,1))<0.001 & ...
        abs(rgbAll(:,2)-rgbUnq(l,2))<0.001 & ...
        abs(rgbAll(:,3)-rgbUnq(l,3))<0.001;
        if subjNum==10
            defocus875predTmp(indStiml) = optDistAll(indStiml)-humanWaveDefocusS10(wvInFocus(l),875);
        end
        if subjNum==1
            defocus875predTmp(indStiml) = optDistAll(indStiml)-humanWaveDefocusS1(wvInFocus(l),875);
        end
        if subjNum==3
            defocus875predTmp(indStiml) = optDistAll(indStiml)-humanWaveDefocusS3(wvInFocus(l),875);
        end
        if subjNum==5
            defocus875predTmp(indStiml) = optDistAll(indStiml)-humanWaveDefocusS5(wvInFocus(l),875);
        end                
        if subjNum==9
            defocus875predTmp(indStiml) = optDistAll(indStiml)-humanWaveDefocusS9(wvInFocus(l),875);
        end
        if subjNum==16
            defocus875predTmp(indStiml) = optDistAll(indStiml)-humanWaveDefocusS16(wvInFocus(l),875);
        end
        if subjNum==17
            defocus875predTmp(indStiml) = optDistAll(indStiml)-humanWaveDefocusS17(wvInFocus(l),875);
        end
        if subjNum==18
            defocus875predTmp(indStiml) = optDistAll(indStiml)-humanWaveDefocusS18(wvInFocus(l),875);
        end                
        if subjNum==20
            defocus875predTmp(indStiml) = optDistAll(indStiml)-humanWaveDefocusS20(wvInFocus(l),875);
        end
        wvInFocusAll(indStiml) = wvInFocus(l);
    end
    
    RMSE = sqrt(mean((defocus875predTmp(:)-defocus875(:)).^2));
else
    optDistUnq = unique(optDistAll);
    defocus875predTmp = zeros([size(rgbUnq,1) length(optDistUnq)]);
    defocus875mean = zeros([size(rgbUnq,1) length(optDistUnq)]);
    for l = 1:size(rgbUnq,1)
        for k = 1:length(optDistUnq)
            indStiml = abs(rgbAll(:,1)-rgbUnq(l,1))<0.001 & ...
                       abs(rgbAll(:,2)-rgbUnq(l,2))<0.001 & ...
                       abs(rgbAll(:,3)-rgbUnq(l,3))<0.001 & ...
                       abs(optDistAll-optDistUnq(k))<0.001;
            defocus875mean(l,k) = mean(defocus875(indStiml));
            if subjNum==10
                defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusS10(wvInFocus(l),875);
            end
            if subjNum==1
                defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusS1(wvInFocus(l),875);
            end
            if subjNum==3
                defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusS3(wvInFocus(l),875);
            end
            if subjNum==5
                defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusS5(wvInFocus(l),875);
            end
            if subjNum==9
                defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusS9(wvInFocus(l),875);
            end
            if subjNum==16
                defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusS16(wvInFocus(l),875);
            end
            if subjNum==17
                defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusS17(wvInFocus(l),875);
            end
            if subjNum==18
                defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusS18(wvInFocus(l),875);
            end
            if subjNum==20
                defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusS20(wvInFocus(l),875);
            end
        end
    end
    % defocus875mean = bsxfun(@minus,defocus875mean,lagLeadTerm);
    RMSE = sqrt(mean((defocus875predTmp(:)-defocus875mean(:)).^2));
end

end
