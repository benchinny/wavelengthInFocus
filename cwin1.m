% present 2 stimuli on 2 screens 
%[RR RG RB] =1.4010    4.1090    0.4064

%[LR LG LB]= 1.3480    3.9330    0.3604

function cwin1(img1, img2, cf, rc, window1, window2);

        tex_1 = Screen('MakeTexture', window1, imcal(circshift(img1, rc(1,:)), cf(:,1)));
        Screen('DrawTexture', window1, tex_1);
        tex_2 = Screen('MakeTexture', window2, imcal(circshift(img2, rc(2,:)), cf(:,2)));
        Screen('DrawTexture', window2, tex_2);
        Screen('Flip', window1);
        Screen('Flip', window2);
        Screen('Close', tex_1); Screen('Close', tex_2);


        %%
