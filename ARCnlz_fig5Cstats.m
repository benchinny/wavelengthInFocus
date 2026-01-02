function pNull = ARCnlz_fig5Cstats(dataPath)

% function for calculating statistic associated with Figure 5C--uses
% bootstrapping to figure out the chances of the peak not being at 0
% (pNull)

rng(1); % SET RANDOM SEED

% LIST OF ALL SUBJECTS
subjNumAll = [1 3 5 10 16 17 18 20];

% INITIALIZE VECTOR FOR STORING PROBABILITY OF PEAK NOT BEING AT 0
pNull = [];

for i = 1:length(subjNumAll) % LOOP OVER SUBJECTS
    % CALL HELPER FUNCTION
    [~,~,~,~,~,~,~,~,~,~,PCboots] = ARCacuityAnalyzeDataOnly(subjNumAll(i),0,dataPath);
    % INITIALIZE VECTOR INDICATING WHETHER OR NOT THE PEAK IS AT 0 FOR EACH
    % RESAMPLING
    bPeakNot0 = [];
    for j = 1:size(PCboots,2) % LOOP OVER RESAMPLES
        % GRAB A RESAMPLING
        PCtmp = PCboots(:,j);
        % DETERMINE WHETHER OR NOT THE MAX IS AT 0
        [PCmaxTmp,~] = max(PCtmp);
        bPeakNot0(j) = PCmaxTmp>PCtmp(5);
    end
    % STORE PROBABILITY
    pNull(i) = 1-mean(bPeakNot0);
end