% %210517 find z0 based on trombone power
function z0=fnz0(a0, ACL)
% LeftTromPwr=a18(5,1);
% RightTromPwr=a18(5,2);


if length(a0)==1
    LeftTromPwr=a0;
    RightTromPwr=a0;
else
    LeftTromPwr=a0(1);
    RightTromPwr=a0(2);
end


 if ACL==0
            %Avg SYSTEM tca values  WITHOUT ACL                                             
            z0= [0	           0	     0	      0     %set default value of z0 as the average values shown here
                -0.91666	2.05	-0.69583	1.3625
                -1.24168	3.3	   -1.10834	   1.44585  ]; 
            
           tc0=Ltca_noACL(LeftTromPwr); %get Left side system tca values WITHOUT ACL based on trombone position
           z0(2,1:2)=tc0(1,:);  %set 
           z0(3,1:2)=tc0(2,:);  %set 
           
           tc0=Rtca_noACL(RightTromPwr); %get Right side system tca values WITHOUT ACL based on trombone position
           z0(2,3:4)=tc0(1,:); %set right 2nd row (Red-green tca) 
           z0(3,3:4)=tc0(2,:) %set right 3rd row (Red-blue tca)
           
%  elseif ACL==1
 else
            %Avg SYSTEM tca values  WITH ACL           
            z0=[ 0	           0	        0	       0
                -0.46667	2.10833 	-0.33332	2.07501
                -0.39999	3.60834	    -0.25834	3.16666 ];
            
           tc0=Ltca_ACL(LeftTromPwr); %get Left side system tca values WITH ACL based on trombone position
           z0(2,1:2)=tc0(1,:);  %set 
           z0(3,1:2)=tc0(2,:);  %set 
           
           tc0=Rtca_ACL(RightTromPwr); %get Right side system tca values WITH ACL based on trombone position
           z0(2,3:4)=tc0(1,:);  %set 
           z0(3,3:4)=tc0(2,:);  %set 
 end


%check the z0 values using power of 8 and 14 on both sides and the values
%obtained match to the data on fitted lines



% % %UPDATE V02: 03/22/21
% % % Link to data:https://docs.google.com/spreadsheets/d/1p_3Ik-r6rub0O3HxF66YCufHOesUDgL5U_p2-u0Cb5c/edit#gid=0
% % %This script calculates system TCA (z0) based on Swati's calibration data
% % %collected on March 17 to 19, 2021. Monochromatic camera with a higher
% % %pixel density was used (along with variable focus lens). It was placed on
% % %x,y,z micrometer stage. Spacer was used between cam n lens. Z location of
% % %camera lens iris dot was 133mm. First cam was adjusted such that image
% % looked centered on the sensor (fov of camera as seen on computer screen),
% % and also looked squared (i.e hor line is horizontal)
% % %Then Camera was always centered laterally before doing 
% % % subjective TCA alignment tasks as seen using the camera, for each trombon
% % % position. Data for 10 different trombone positions was collected.
% % %For each position, total of 12 alignment tasks were done.
% %z location of camera is decided based on where the image on camera disappears
% %uniformly on laterally moving the camera
% % 
% % % Steps: Set the trombone position, center the camera, do the tca
% % alignment (repeat on each side, with and without ACL)
% % 
% % %Initially, set the system tca, z0 as the average value of the data
% % %obtained (instead of setting z0 to zeros)
% % 
% % %There are 4 functions which calculate system tca for each trombone power/position,
% % %based on line fit equations to the data collected in different conditions
% % %Rtca_noACL(RightTromPwr)  Rtca_ACL(RightTromPwr) Ltca_noACL(LeftTromPwr) Ltca_ACL(LeftTromPwr)
% % %the above functions use the data fitting shown in the following files
% % %RightNOACL.m   RightWITHACL.m   LeftNOACL.m   LeftWITHACL.m
% % %Format of z0 is: 
% % %z0= [Red(leftY leftX rightY rightX);
% % %    Green(leftY leftX rightY rightX);
% % %    Blue(leftY leftX rightY rightX)]
% 
% %JUST LOAD LCA saved data file to get trombone power
% LeftTromPwr=a18(5,1);
% RightTromPwr=a18(5,2);
% 
% 
%  if ACL==0
%             %Avg SYSTEM tca values  WITHOUT ACL                                             
%             z0= [0	           0	     0	      0     %set default value of z0 as the average values shown here
%                 -0.91666	2.05	-0.69583	1.3625
%                 -1.24168	3.3	   -1.10834	   1.44585  ]; 
%             
%            tc0=Ltca_noACL(LeftTromPwr); %get Left side system tca values WITHOUT ACL based on trombone position
%            z0(2,1:2)=tc0(1,:);  %set 
%            z0(3,1:2)=tc0(2,:);  %set 
%            
%            tc0=Rtca_noACL(RightTromPwr); %get Right side system tca values WITHOUT ACL based on trombone position
%            z0(2,3:4)=tc0(1,:); %set right 2nd row (Red-green tca) 
%            z0(3,3:4)=tc0(2,:) %set right 3rd row (Red-blue tca)
%            
%  elseif ACL==1
%             %Avg SYSTEM tca values  WITH ACL           
%             z0=[ 0	           0	        0	       0
%                 -0.46667	2.10833 	-0.33332	2.07501
%                 -0.39999	3.60834	    -0.25834	3.16666 ];
%             
%            tc0=Ltca_ACL(LeftTromPwr); %get Left side system tca values WITH ACL based on trombone position
%            z0(2,1:2)=tc0(1,:);  %set 
%            z0(3,1:2)=tc0(2,:);  %set 
%            
%            tc0=Rtca_ACL(RightTromPwr); %get Right side system tca values WITH ACL based on trombone position
%            z0(2,3:4)=tc0(1,:);  %set 
%            z0(3,3:4)=tc0(2,:)  %set 
%  end
% 
% 
% %check the z0 values using power of 8 and 14 on both sides and the values
% %obtained match to the data on fitted lines



%UPDATE V01: 03/22/21
%      if strcmp(ey,'Right')
%          
%         if ACL==0 
%                                                      %z0 is the system TCA
%              z0= [0	           0	     0	     0   %set default value of z0 as the average values shown here
%                 -0.91666	2.05	-0.69583	1.3625
%                 -1.24168	3.3	   -1.10834	   1.44585  ]; %average values from march17-19 calibration with monochromatic camera which had higher pixel density)
%               
%                                                               
%              tc0=Rtca_noACL(RightTromPwr); %get right side system tca values (NO ACL) based on trombone position
%              z0(2,3:4)=tc0(1,:); %update right 2nd row (Red-green tca) 
%              z0(3,3:4)=tc0(2,:); %update right 3rd row (Red-blue tca)
%          
%         elseif ACL==1
%            
%              
%             z0=[ 0	           0	        0	       0     %Avg values with WITH ACL
%                 -0.46667	2.10833 	-0.33332	2.07501
%                 -0.39999	3.60834	    -0.25834	3.16666 ]; %average values from march17-19 calibration with monochromatic camera which had higher pixel density)
%               
%             tc0=Rtca_ACL(RightTromPwr); %get right side system tca values (WITH ACL) based on trombone position
%             z0(2,3:4)=tc0(1,:);  %set 
%             z0(3,3:4)=tc0(2,:);  %set 
%         end
%         
%       elseif strcmp(ey,'Left')
%          
%           if ACL==0
%               
%             z0= [0	           0	     0	      0
%                 -0.91666	2.05	-0.69583	1.3625
%                 -1.24168	3.3	   -1.10834	   1.44585  ]; %average values from march17-19 calibration with monochromatic camera which had higher pixel density)
%            tc0=Ltca_noACL(LeftTromPwr); %get Left side system tca values (WITHOUT ACL) based on trombone position
%            z0(2,1:2)=tc0(1,:);  %set 
%            z0(3,1:2)=tc0(2,:);  %set 
%            
%           elseif ACL==1
%             
%             z0=[ 0	           0	        0	       0
%                 -0.46667	2.10833 	-0.33332	2.07501
%                 -0.39999	3.60834	    -0.25834	3.16666 ]; %average values from march17-19 calibration with monochromatic camera which had higher pixel density)
%            tc0=Ltca_ACL(LeftTromPwr); %get Left side system tca values (WITH ACL) based on trombone position
%            z0(2,1:2)=tc0(1,:);  %set 
%            z0(3,1:2)=tc0(2,:);  %set 
%           end
%           
%      elseif strcmp(ey,'Stereo')
%          
%           if ACL==0
%               
%             z0= [0	           0	     0	      0
%                 -0.91666	2.05	-0.69583	1.3625
%                 -1.24168	3.3	   -1.10834	   1.44585  ]; %average values from march17-19 calibration with monochromatic camera which had higher pixel density)
%            tc0=Ltca_noACL(LeftTromPwr); %get Left side system tca values (WITHOUT ACL) based on trombone position
%            z0(2,1:2)=tc0(1,:);  %set 
%            z0(3,1:2)=tc0(2,:);  %set 
%             tc0=Rtca_noACL(RightTromPwr); %get right side system tca values (NO ACL) based on trombone position
%              z0(2,3:4)=tc0(1,:); %update right 2nd row (Red-green tca) 
%              z0(3,3:4)=tc0(2,:); %update right 3rd row (Red-blue tca)
%            
%           elseif ACL==1
%             
%             z0=[ 0	           0	        0	       0
%                 -0.46667	2.10833 	-0.33332	2.07501
%                 -0.39999	3.60834	    -0.25834	3.16666 ]; %average values from march17-19 calibration with monochromatic camera which had higher pixel density)
%            tc0=Ltca_ACL(LeftTromPwr); %get Left side system tca values (WITH ACL) based on trombone position
%            z0(2,1:2)=tc0(1,:);  %set 
%            z0(3,1:2)=tc0(2,:);  %set 
%            tc0=Rtca_ACL(RightTromPwr); %get right side system tca values (WITH ACL) based on trombone position
%             z0(2,3:4)=tc0(1,:);  %set 
%             z0(3,3:4)=tc0(2,:);  %set 
%           end
%           
%      end