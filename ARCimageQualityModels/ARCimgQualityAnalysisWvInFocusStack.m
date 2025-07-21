%%

wvInFocusSTall = [];
wvInFocusXCall = [];
for i = 1:4
    [wvInFocusST, wvInFocusXC] = ARCimgQualityAnalysisVarPilot2abb(i);
    wvInFocusSTall(:,i) = wvInFocusST;
    wvInFocusXCall(:,i) = wvInFocusXC;
end

save('/Users/benjaminchin/Documents/ARchromaScraps/wvInFocusComparison.mat','wvInFocusSTall','wvInFocusXCall');