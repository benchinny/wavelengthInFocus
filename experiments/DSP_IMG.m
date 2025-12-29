       




%%211102 DSP_IMG

clear LCAim
[window1, window2, vbl0]=strt_psych0(screenNumber-2, screenNumber-1, 0);

        %%input a output b
        cf=ones(3,2);

        %initial  power=13.5 AR    
        %LCAim='texture0_1080_newfill_malt.png';    
        %im2R='black.png';       
         im2L='texture0_1080_newfill_malt.png';
        im2R=im2L ;
        
            deg=-3;            
            zaber(name_map('rotation')).move_deg(deg); %%-6400            
            zaber(name_map('rotation')).control.getposition  
            
            if exist('sr') ~=  1; sr=[0 0]; end
%         dmnd=[-0.5:0.5:3]
           opto(name_map('r_disp')).control.setFocalPower(14.4+sr(2));% -dmnd(k0));
           opto(name_map('r_disp')).control.getFocalPower.focal_power
           
           opto(name_map('l_disp')).control.setFocalPower(14+sr(1));%-dmnd(k0));
           opto(name_map('l_disp')).control.getFocalPower.focal_power

        [iLf iRf]=cwin3(imread(im2L), imread(im2R) , cf, rc00, window2, window1);
%         for k0=1:length(dmnd)
  %            disp(n2s(dmnd(k0)))
%         KbWait([], 2);
%         end
    %img=imread(LCAim); %% image to show
    %wn=cwin0(img, 'Stereo', cf, rc00, window1, window2); 
    
    
    
        KbWait([], 2);
        [iLf iRf]=cwin3(imread("black.png"), imread("black.png") , cf, rc00, window2, window1);

        
        clear LCAim
        sca   







