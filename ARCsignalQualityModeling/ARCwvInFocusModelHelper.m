function [RMSE, defocus875mean, defocus875predTmp, rgbUnq, optDistUnq] = ARCwvInFocusModelHelper(subjNum,defocus875,rgbAll,optDistAll,w,sigQualType,dataPath)

% LIST OF ALL SUBJECTS
subjNumListAll = [1 3 5 10 16 17 18 20];
% FIND subjNum POSITION IN ARRAY
indLCA = find(subjNumListAll==subjNum);
% LOAD PRE-SAVED LCA PARAMETERS
load(fullfile(dataPath,'data','PresavedFigureData','LCAparams.mat'),'q1bestAll','q2bestAll','q3bestAll');
q1 = q1bestAll(indLCA);
q2 = q2bestAll(indLCA);
q3 = q3bestAll(indLCA);

rgbUnq = unique(rgbAll,'rows'); % UNIQUE COLOR CONDITIONS
wvInFocus = zeros([size(rgbUnq,1) 1]); % INITIALIZE WAVELENGTH IN FOCUS VECTOR
for l = 1:size(rgbUnq,1) % FOR EACH COLOR CONDITION
    % GET THE WAVELENGTH THAT SHOULD BE IN FOCUS FOR EACH STIMULUS
    if strcmp(sigQualType,'xcorr')
       wvInFocus(l,:) = ARCwvInFocusConesMeanZspatFilter(subjNum,l,w,dataPath);
    elseif strcmp(sigQualType,'strehl')
       wvInFocus(l,:) = ARCwvInFocusConesMeanZstrehl(subjNum,rgbUnq(l,:),w,dataPath);
    elseif strcmp(sigQualType,'deltapass')
       wvInFocus(l,:) = ARCwvInFocusConesMeanZdeltaPass(subjNum,rgbUnq(l,:),w,dataPath);
    else
        error('ARCwvInFocusModelHelper: invalid value of parameter sigQualType. Must be xcorr or strehl.');
    end
    display(['stim ' num2str(l)]);
end

optDistUnq = unique(optDistAll); % UNIQUE STIM DISTANCES
% INITIALIZE MATRICES FOR STORING PREDICTIONS AND SORTING MEASUREMENTS
defocus875predTmp = zeros([size(rgbUnq,1) length(optDistUnq)]);
defocus875mean = zeros([size(rgbUnq,1) length(optDistUnq)]);
for l = 1:size(rgbUnq,1) % LOOP OVER COLOR CONDITIONS
    for k = 1:length(optDistUnq) % LOOP OVER STIM DISTANCES
        indStiml = abs(rgbAll(:,1)-rgbUnq(l,1))<0.001 & ...
                   abs(rgbAll(:,2)-rgbUnq(l,2))<0.001 & ...
                   abs(rgbAll(:,3)-rgbUnq(l,3))<0.001 & ...
                   abs(optDistAll-optDistUnq(k))<0.001;
        % GET MEAN DEFOCUS ABERRATION IN A CONDITION AFTER SORTING
        defocus875mean(l,k) = mean(defocus875(indStiml));
        defocus875predTmp(l,k) = optDistUnq(k)-humanWaveDefocusParameterizedARC(wvInFocus(l),875,q1,q2,q3);
    end
end

RMSE = sqrt(mean((defocus875predTmp(:)-defocus875mean(:)).^2));

end
