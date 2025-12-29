      
%%211102 DSP_IMG

clear LCAim
[window1, window2, vbl0]=strt_psych0(screenNumber-1, screenNumber, 0);

        %%input a output b
        cf=ones(3,2);

        %initial  power=13. 5 AR
        %LCAim='texture0_1080_newfill_malt.png';    
        im2R='black.png';       
        % im2L='texture0_1080_newfill_malt.png';
        im2L=im2R;
         
        
    %img=imread(LCAim); %% image to show
    %wn=cwin0(img, 'Stereo', cf, rc00, window1, window2); 
    [iLf iRf]=cwin3(imread(im2L), imread(im2R) , cf, rc00, window2, window1);
        KbWait([], 2);
        clear LCAim
        sca








         