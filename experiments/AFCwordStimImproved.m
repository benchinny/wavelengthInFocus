function imgWord = AFCwordStimImproved(word2display,imSzPix,wordColor)

% function imgWord = AFCwordStimImproved(word2display,imSize)
%
% makes an image of a word
%
% inputs: 
%         word2display: string of characters forming a word [1 x n]
%         imSzPix     : image size in pixels [1 x 2]
%         charLocations: locations of characters [n x 2]

   imgWord = zeros(imSzPix);
   imgWord = insertText(imgWord,imSzPix./2,word2display,'TextColor',wordColor,'BoxColor','black','FontSize',round((imSzPix(1)*0.5)),'AnchorPoint','Center'); 
   imgWord = imgWord.*255;
end