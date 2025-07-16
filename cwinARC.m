function []=cwinARC(img1, img2, window1, window2)
%         iAlign = zeros([960 960 4]);
%         iAlign(:,:,1) = 100;
%         iAlign(430:530,281:664,2) = 100;
%         iAlign(:,:,4) = 0.3;
        tex_1 = Screen('MakeTexture', window1, uint8(img1));
        Screen('DrawTexture', window1, tex_1);
        tex_2 = Screen('MakeTexture', window2, uint8(img2));
        Screen('DrawTexture', window2, tex_2);
%        texAlign = Screen('MakeTexture', window1, iAlign);
%        Screen('DrawTexture', window1, texAlign);
        Screen('Flip', window1);
        Screen('Flip', window2);
        Screen('Close', tex_1); Screen('Close', tex_2); 
%        Screen('Close', texAlign);
