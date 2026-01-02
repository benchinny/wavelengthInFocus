function [meanLag, CI95lag, LtoMratioMean, LtoMratioCI95, StoLplusMratioMean, StoLplusMratioCI95] = ARCnlz_Fig4statsHelper(pFitAll,wLpropMinAll,wLMminAll,wS)

% helper function for calculating parameters stats reported in the
% manuscript. 
% 
% pFitAll: presaved lag/lead parameters
% wLpropMinAll: L/M cone weight ratio
% wLMminAll: weight on (L+M) during fit (does not include weight on S
%            initially)
%
% meanLag: mean lag value
% CI95lag: 95% CIs for lags
% LtoMratioMean: mean L to M weight ratio
% LtoMratioCI95: 95% CI for L to M weight ratio
% StoLplusMratioMean: mean S/(L+M) ratio
% StoLplusMratioCI95: 95% CIs on (L+M)-S ratio

% LIST OF SUBJECTS
subjNumAll = [1 3 5 10 16 17 18 20];
% INITIALIZE ARRAYS FOR STORING (L+M)/S WEIGHT RATIO (IN CASE NOT USING)
StoLplusMratioMean = []; 
StoLplusMratioCI95 = [];

% CALCULATE MEAN AND 95% CIS FOR L/M WEIGHT RATIO
LtoMratioMean = mean(wLpropMinAll);
LtoMratioCI95 = 1.96.*std(wLpropMinAll)./sqrt(length(subjNumAll));
% CALCULATE MEAN AND 95% CIS FOR (L+M)/S WEIGHT RATIO
if ~isempty(wLMminAll) && ~isempty(wS)
    StoLplusMratioMean = mean(wS./wLMminAll);
    StoLplusMratioCI95 = 1.96.*std(wS./wLMminAll)./sqrt(length(subjNumAll));
end

% LIST OF OPTICAL DISTANCES USED IN EXPERIMENT
optDist = [1.5 2.5 3.5];

% ARRAY REPRESENTING HOW MUCH THE FITS NEED TO BE SHIFTED FOR EACH SUBJECT
% AND DISTANCE
defocusAt875fit = [];

for i = 1:size(pFitAll,2) % LOOP OVER SUBJECTS
    for j = 1:length(optDist) % LOOP OVER STIMULUS DISTANCE
        % CALCULATE LAG BASED ON LAG/LEAD PARAMETERS FROM FIT
        defocusAt875fit(i,j) = -optDist(j)*pFitAll(1,i) - pFitAll(2,i);
    end
end

% CALCULATE MEAN AND 95% CIS FOR LAGS
meanLag = mean(defocusAt875fit);
CI95lag = 1.96.*std(defocusAt875fit,0,1)./sqrt(8);

end