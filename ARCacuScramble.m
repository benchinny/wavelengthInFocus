%%

rgb = [0.569 0 1];
meanFocstmOptDst = [2.5]*1.2255;
focStmOptDstIncr = -1.2:0.3:1.2;
trlPerLvl = 6;

rgbAll = [];
meanFocstmOptDstAll = [];
focStmOptDstIncrAll = [];
indAcuRGBall = [];
indScramble = [];
maskBrightness = 0;
maskSize = [100 100];
gammaR = 2.5;
gammaG = 2.7;
gammaB = 2.3;
frqCpd = 15;
rgbAcuRGB = [0.569 0 0; 0 0.432 0; 0 0 1];

% CAREFUL ATTEMPT TO BLOCK CONDITIONS SO EACH OPTICAL DISTANCE INCREMENT IS
% PRESENTED ONCE PER BLOCK
for i = 1:size(rgb,1)
   for j = 1:length(meanFocstmOptDst)
       for m = 1:size(rgbAcuRGB,1)
           for k = 1:length(focStmOptDstIncr)           
               rgbAll(end+1,:) = rgb(i,:);
               meanFocstmOptDstAll(end+1,:) = meanFocstmOptDst(j);
               focStmOptDstIncrAll(end+1,:) = focStmOptDstIncr(k);
               indAcuRGBall(end+1,:) = m;
           end
       end
       for l = 1:trlPerLvl
          indScramble = [indScramble; randperm(length(focStmOptDstIncr)*size(rgbAcuRGB,1))'];
       end
   end
end

% RANDOMIZING TRIALS
indScramble = indScramble+imresize(length(focStmOptDstIncr).*size(rgbAcuRGB,1).*[0:(trlPerLvl*size(rgb,1)*length(meanFocstmOptDst)-1)]',size(indScramble),'nearest');
rgbAll = repmat(rgbAll,[trlPerLvl 1]);
meanFocstmOptDstAll = repmat(meanFocstmOptDstAll,[trlPerLvl 1]);
focStmOptDstIncrAll = repmat(focStmOptDstIncrAll,[trlPerLvl 1]);
indAcuRGBall = repmat(indAcuRGBall,[trlPerLvl 1]);
%%
rgbAll = rgbAll(indScramble,:);
meanFocstmOptDstAll = meanFocstmOptDstAll(indScramble);
focStmOptDstIncrAll = focStmOptDstIncrAll(indScramble);
indAcuRGBall = indAcuRGBall(indScramble);

% ADD DUMMY TRIAL RIGHT AT THE END (PECULIAR TO WAY CODE IS WRITTEN)
rgbAll(end+1,:) = [0 0 0];
focStmOptDstIncrAll(end+1,:) = 0;
meanFocstmOptDstAll(end+1,:) = 3;
stimSizePixAll(end+1,:) = 10;
indAcuRGBall(end+1,:) = 1;
