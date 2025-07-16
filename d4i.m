%% put d0i into a forloop


function [i4 i3 i2 i1 i0]=d4i(x, ms);
rc=size(x);
i0=zeros(rc);
i1=i0; i2=i0; i3=i0; i4=i0;

for r=1:rc(1);
    for c=1:rc(2);
        [y4 y3 y2 y1 y0]=d0i(x(r,c), ms); 
        i4(r,c)=y4;
        i3(r,c)=y3;
        i2(r,c)=y2;
        i1(r,c)=y1;
        i0(r,c)=y0;
        
    end
end