function [y5 y4 y3 y1 y0]=d0i(x, ms);
% number to integer and decimal
% based on d2i, including correction for negative values.
% [negative=1 correct_decimal decimal integer number]
% x=-5.7
% ms=10;

y0=round(ms.*x)./ms; % correct number to get exact;    
y1=fix(y0); %integer
y2=round(ms.*y0)-y1.*ms; % get first decimal
y3=round(abs(y2)); % decimal point
v0=[1 10 9 8 7 6 5 4 3 2];
if y0>=0;
    y4=y3+1;
    y5=0;
else
   y4=v0(y3+1);
   y5=1;
end