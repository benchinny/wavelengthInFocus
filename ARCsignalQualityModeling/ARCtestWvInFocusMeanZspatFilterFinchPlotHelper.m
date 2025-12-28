function [RMSE, defocus550, defocus550predTmp, rbRatio, optDistUnq] = ARCtestWvInFocusMeanZspatFilterFinchPlotHelper(subjNum,defocus550,rbRatio,optDistUnq,w)

wvInFocus = zeros([length(rbRatio) 6]);
for l = 5
    wvInFocus(l,:) = ARCwvInFocusConesMeanZspatFilterFinch(subjNum,l,w);
    display(['stim ' num2str(l)]);
end 
indBad = wvInFocus<400 | wvInFocus>720;
wvInFocus(indBad) = NaN;

defocus550predTmp = zeros([length(rbRatio) 6]);
for l = 1:size(rbRatio,2)
    if subjNum==10
        defocus550predTmp = optDistUnq-humanWaveDefocusS10(wvInFocus,550);
    end
    if subjNum==1
        defocus550predTmp = optDistUnq-humanWaveDefocusS1(wvInFocus,550);
    end
    if subjNum==3
        defocus550predTmp = optDistUnq-humanWaveDefocusS3(wvInFocus,550);
    end
    if subjNum==5
        defocus550predTmp = optDistUnq-humanWaveDefocusS5(wvInFocus,550);
    end
    if subjNum==9
        defocus550predTmp = optDistUnq-humanWaveDefocusS9(wvInFocus,550);
    end
    if subjNum==16
        defocus550predTmp = optDistUnq-humanWaveDefocusS16(wvInFocus,550);
    end
    if subjNum==17
        defocus550predTmp = optDistUnq-humanWaveDefocusS17(wvInFocus,550);
    end
    if subjNum==18
        defocus550predTmp = optDistUnq-humanWaveDefocusS18(wvInFocus,550);
    end
    if subjNum==20
        defocus550predTmp = optDistUnq-humanWaveDefocusS20(wvInFocus,550);
    end
    if subjNum==21
        defocus550predTmp = optDistUnq-humanWaveDefocusS21(wvInFocus,550);
    end
end
RMSE = sqrt(mean((defocus550predTmp(:)-defocus550(:)).^2));

end
