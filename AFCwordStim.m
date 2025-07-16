function imgWord = AFCwordStim(word2display,imSzPix,charLocations,wordColor,fontSize)

% function imgWord = AFCwordStim(word2display,imSize)
%
% makes an image of a word
%
% inputs: 
%         word2display: string of characters forming a word [1 x n]
%         imSzPix     : image size in pixels [1 x 2]
%         charLocations: locations of characters [n x 2]

if length(word2display)~=size(charLocations,1)
    error('AFCwordStim: length of word2display must be equal to number of rows in charLocations!');
end

   imgWord = zeros(imSzPix);
   for i = 1:length(word2display)
       imgWord = insertText(imgWord,charLocations(i,:),word2display(i),'TextColor',wordColor,'BoxColor','black','FontSize',fontSize); 
   end
   imgWord = imgWord.*255;
end