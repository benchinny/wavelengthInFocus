function [im2L0, im2L1, im2R0, im2R1] = AFCtcaStmImg(e0, e1, bxyL, bxyR)

global sz cf rc00

ms=10;
[j4L j3L j2L j1L j0L]=d4i(bxyL, ms); % system+subject TCA
[j4R j3R j2R j1R j0R]=d4i(bxyR, ms); % system+subject TCA

i1=ones(sz); i255=255.*i1;

irL0= embd(circshift(e0{j3L(1,1), j3L(1,2)}    ,[j1L(1,1)-j4L(1,1) j1L(1,2)-j4L(1,2)] ),i255); %window2 LFT screen RED

igL0= embd(circshift(e0{j3L(2,1), j3L(2,2)}    ,[j1L(2,1)-j4L(2,1) j1L(2,2)-j4L(2,2)] ),i255); %window2 LFT screen GREEN

ibL0= embd(circshift(e0{j3L(3,1), j3L(3,2)}    ,[j1L(3,1)-j4L(3,1) j1L(3,2)-j4L(3,2)] ),i255); %window2 LFT screen BLUE
     
irR0= embd(circshift(e0{j3R(1,1), j3R(1,2)}    ,[j1R(1,1)-j4R(1,1) j1R(1,2)-j4R(1,2)] ),i255); %window2 RGT screen RED

igR0= embd(circshift(e0{j3R(2,1), j3R(2,2)}    ,[j1R(2,1)-j4R(2,1) j1R(2,2)-j4R(2,2)] ),i255); %window2 RGT screen GREEN

ibR0= embd(circshift(e0{j3R(3,1), j3R(3,2)}    ,[j1R(3,1)-j4R(3,1) j1R(3,2)-j4R(3,2)] ),i255); %window2 RGT screen BLUE
     
irL1= embd(circshift(e1{j3L(1,1), j3L(1,2)}    ,[j1L(1,1)-j4L(1,1) j1L(1,2)-j4L(1,2)] ),i255); %window2 LFT screen RED

igL1= embd(circshift(e1{j3L(2,1), j3L(2,2)}    ,[j1L(2,1)-j4L(2,1) j1L(2,2)-j4L(2,2)] ),i255); %window2 LFT screen GREEN

ibL1= embd(circshift(e1{j3L(3,1), j3L(3,2)}    ,[j1L(3,1)-j4L(3,1) j1L(3,2)-j4L(3,2)] ),i255); %window2 LFT screen BLUE 
 
irR1= embd(circshift(e1{j3R(1,1), j3R(1,2)}    ,[j1R(1,1)-j4R(1,1) j1R(1,2)-j4R(1,2)] ),i255); %window2 RGT screen RED

igR1= embd(circshift(e1{j3R(2,1), j3R(2,2)}    ,[j1R(2,1)-j4R(2,1) j1R(2,2)-j4R(2,2)] ),i255); %window2 RGT screen GREEN

ibR1= embd(circshift(e1{j3R(3,1), j3R(3,2)}    ,[j1R(3,1)-j4R(3,1) j1R(3,2)-j4R(3,2)] ),i255); %window2 RGT screen BLUE

im2L0=cat(3, irL0, igL0, ibL0); 
im2R0=cat(3, irR0, igR0, ibR0); 
          
im2L1=cat(3, irL1, igL1, ibL1); 
im2R1=cat(3, irR1, igR1, ibR1); 

end