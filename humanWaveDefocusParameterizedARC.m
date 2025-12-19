function D = humanWaveDefocusParameterizedARC(waveRef,wave,q1,q2,q3)
% Defocus in diopters as a function of wavelength
%
% Syntax:
%   D = humanWaveDefocusParameterizedARC(waveRef,wave,q1,q2,q3)
%
% Description:
%    This is a function fit to the data from Bedford and Wyszecki and Wald
%    on human chromatic aberration.
%
% Inputs:
%    wave - Vector. Wavelength vector, in nanometers.
%
% Outputs:
%    D    - Vector. Defocus vector.
%
% Optional key/value pairs:
%    None.

% Constants for formula to compute defocus in diopters (D) as a function of
% wavelength for human eye. From 
% Larry N. Thibos, Ming Ye, Xiaoxiao Zhang, and Arthur Bradley, 
% "The chromatic eye: a new reduced-eye model of ocular chromatic aberration 
% in humans," Appl. Opt. 31, 3594-3600 (1992)

% This is the human defocus as a function of wavelength. This formula
% converts the wave in nanometers to wave in microns. D is in diopters.
Dref = q1 - (q2 ./ (waveRef * 1e-3 - q3));
D = q1 - (q2 ./ (wave * 1e-3 - q3))-Dref;

return
