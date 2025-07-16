function OTF = ARCfreqFilterGauss2D(nSmpXY,szDegXY,muFreq,stdFreq)

% SAMPLES PER DEGREE
smpPerDegXY = nSmpXY./szDegXY;

% NYQUIST FREQUENCY
nyqFrqXY = smpPerDegXY./2;

% OBTAIN FREQUENCIES
frqX = linspace(-nyqFrqXY(1),nyqFrqXY(1),nSmpXY(1));
frqY = linspace(-nyqFrqXY(2),nyqFrqXY(2),nSmpXY(2));

% CONVERT TO MESHGRID
[frqXX, frqYY] = meshgrid(frqX,frqY);

OTF = normpdf(sqrt(frqXX.^2+frqYY.^2),muFreq,stdFreq);
OTF = OTF./max(OTF(:));

end