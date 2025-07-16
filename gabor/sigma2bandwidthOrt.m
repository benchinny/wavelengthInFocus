function [BWort] = sigma2bandwidthOrt(frqCpd, sigmaDeg)
  
% function [sigmaDeg sigmaCpd] = sigma2bandwidthOrt(frqCpd, sigmaDeg)
%
%   example call: sigmaDeg = sigma2bandwidthOrt(4, .2)
%
% see derivation in Proof_GaussianBandwidth*.doc in ../VisionNotes/
% 
% returns standard deviation in low pass direction of 2D gabor 
% given an orientation bandwidth in radians
%
% frqCpd:   carrier frequency in cycles/deg
% sigmaDeg: standard deviation of Gaussian envelope 
%           in low pass direction of space domain (deg) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BWort:    orientation (linear) bandwidth is the angle (radians) subtended
%           by the width of the gaussian envelope in the low pass direction
%           of frequency space
%
% ***                   see sigma2bandwidthOrt.m                       ***

BWort = 2.*atan(sqrt(log(4))./(2.*pi.*sigmaDeg.*frqCpd));


