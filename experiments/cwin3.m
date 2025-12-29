% present 2 stimuli on 2 screens 
%[RR RG RB] =1.4010    4.1090    0.4064

%[LR LG LB]= 1.3480    3.9330    0.3604

function [i1 i2]=cwin3(img1, img2, cf, rc, window1, window2);
        i1=imcal(circshift(img1, rc(1,:)), cf(:,1));
%         iAlign = zeros([960 960 4]);
%         iAlign(:,:,1) = 100;
%         iAlign(430:530,281:664,2) = 100;
%         iAlign(:,:,4) = 0.3;
        tex_1 = Screen('MakeTexture', window1, i1);
        Screen('DrawTexture', window1, tex_1);
        i2= imcal(circshift(img2, rc(2,:)), cf(:,2));
        tex_2 = Screen('MakeTexture', window2, i2);
        Screen('DrawTexture', window2, tex_2);
%        texAlign = Screen('MakeTexture', window1, iAlign);
%        Screen('DrawTexture', window1, texAlign);
        Screen('Flip', window1);
        Screen('Flip', window2);
        Screen('Close', tex_1); Screen('Close', tex_2); 
%        Screen('Close', texAlign);


        %%
