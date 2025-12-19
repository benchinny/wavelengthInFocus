function wave = humanWaveDefocusInvertParameterizedARC(waveRef,D,q1,q2,q3)
% wavelength as a function of diopters
%
% Syntax:
%   wave = humanWaveDefocusInvertParameterizedARC(waveRef,D,q1,q2,q3)
%
% Description:
%    This is an inversion of a function fit to the data from 
%    Bedford and Wyszecki and Wald on human chromatic aberration.
%
% Inputs:
%    D    - Vector. Defocus vector.
%
% Outputs:
%    wave - Vector. Wavelength vector, in nanometers.
%
% Optional key/value pairs:
%    None.
%

% Constants for formula to compute defocus in diopters (D) as a function of
% wavelength for human eye. From 
% Larry N. Thibos, Ming Ye, Xiaoxiao Zhang, and Arthur Bradley, 
% "The chromatic eye: a new reduced-eye model of ocular chromatic aberration 
% in humans," Appl. Opt. 31, 3594-3600 (1992)

% This is the wavelength as a function of defocus. 
Dref = q1 - (q2 ./ (waveRef * 1e-3 - q3));
q1 = q1-Dref;
wave = ((q2./(q1-D))+q3)./0.001;

return
