function BWoct = sigma2bandwidthOct(frqCpd, sigmaDeg)

% function BWoct = sigma2bandwidthOct(frqCpd, sigmaDeg)
%
%   example call: BWoct = sigma2bandwidthOct(2, 1)
%
% octave bandwidth of a gabor function given the carrier frequency 
% and the standard deviation of the gaussian envelope in space
%
% see derivation in Proof_GaussianBandwidth_*.doc in ../VisionNotes/
%
% frqCpd:    gabor carrier frequency                          (cyc/deg)
% sigmaDeg: standard deviation of spatial gaussian envelope  (deg)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BWoct:    octave bandwidth is log2( fHi / fLo ) 
%           of the gabor's amplitude spectrum
%
% ***                   see bandwidthOct2sigma.m                       ***

% CONVERT SD IN SPACE TO SD IN FREQUENCY
sigmaCpd = 1./(2.*pi.*sigmaDeg);

% HI AND LO FREQUENCYS AT FWHM
fHiCpd = frqCpd + sigmaCpd.*sqrt(log(4));
fLoCpd = frqCpd - sigmaCpd.*sqrt(log(4));

% fHiCpd = frqCpd + sqrt(log(4))./(2.*pi.*sigmaDeg);
% fLoCpd = frqCpd - sqrt(log(4))./(2.*pi.*sigmaDeg);

% OCTAVE BANDWIDTH
BWoct = log2( fHiCpd ./ fLoCpd );

% BWoct = log2( ( frqCpd + sqrt(log(4))./(2.*pi.*sigmaDeg)) ./ ( frqCpd - sqrt(log(4))./(2.*pi.*sigmaDeg)) );