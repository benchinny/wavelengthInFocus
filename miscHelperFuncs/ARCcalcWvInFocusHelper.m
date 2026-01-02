function [defocus875mean, rgbUnq, optDistUnq] = ARCcalcWvInFocusHelper(defocus875,rgbAll,optDistAll)

% helper function for calculating wavelength in focus. Basically takes in 3
% trial-by-trials vectors of defocus values, rgb values, and stimulus
% distances. 

rgbUnq = unique(rgbAll,'rows'); % UNIQUE COLOR CONDITIONS
optDistUnq = unique(optDistAll); % UNIQUE STIM DISTANCES

% INITIALIZE MATRICES FOR STORING PREDICTIONS AND SORTING MEASUREMENTS
defocus875mean = zeros([size(rgbUnq,1) length(optDistUnq)]);
for l = 1:size(rgbUnq,1) % LOOP OVER COLOR CONDITIONS
    for k = 1:length(optDistUnq) % LOOP OVER STIM DISTANCES
        indStiml = abs(rgbAll(:,1)-rgbUnq(l,1))<0.001 & ...
                   abs(rgbAll(:,2)-rgbUnq(l,2))<0.001 & ...
                   abs(rgbAll(:,3)-rgbUnq(l,3))<0.001 & ...
                   abs(optDistAll-optDistUnq(k))<0.001;
        % GET MEAN DEFOCUS ABERRATION IN A CONDITION AFTER SORTING
        defocus875mean(l,k) = mean(defocus875(indStiml));
    end
end

end
