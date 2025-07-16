function [sigmaDeg, sigmaCpd] = bandwidthOrt2sigma(f0cpd, BWort)
  
% function [sigmaDeg sigmaCpd] = bandwidthOrt2sigma(f0cpd, BWort)
%
%   example call: sigmaDeg = bandwidthOrt2sigma(4, 42.*pi./180)
%
% see derivation in Proof_GaussianBandwidth*.doc in ../VisionNotes/
% 
% returns standard deviation in low pass direction of 2D gabor 
% given an orientation bandwidth in radians
%
% f0cpd:    carrier frequency in cycles/deg
% BWort:    orientation bandwidth in radians 
%           the arc subtended by the width of the gaussian envelope 
%           in the low pass direction of frequency space
% %%%%%%%%%%%%%%%%%%%%%%
% sigmaDeg: gaussian standard deviation in the SPACE     domain (deg) 
% sigmaCpd: gaussian standard deviation in the FREQUENCY domain (cpd) 
%
% ***                   see sigma2bandwidthOct.m                       ***

sigmaDeg = sqrt(log(4))./(2.*pi.*f0cpd.*tan(BWort./2));
sigmaCpd = 1./(2.*pi.*sigmaDeg);

