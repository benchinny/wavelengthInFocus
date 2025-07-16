function RMSE = ARCnlzSwitchingFunc(deltaA,deltaRB1,deltaRB2,deltaS,p)

d=p(1);
wS=p(2);
rbThreshold=p(3);

c = zeros(size(deltaS));
wS = wS.*ones(size(deltaS));

% COMPARE DIFFERENCES BETWEEN RED AND BLUE LUMINANCES TO THRESHOLD
% THEN ADD A POSITIVE OR NEGATIVE CONSTANT
c(deltaRB1<rbThreshold & deltaRB2>rbThreshold)=d;
c(deltaRB1>rbThreshold & deltaRB2<rbThreshold)=-d;

deltaApredicted = wS.*deltaS + c;
RMSE = sqrt(mean((deltaApredicted-deltaA).^2));

end