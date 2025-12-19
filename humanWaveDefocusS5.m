function D = humanWaveDefocusS5(waveRef,wave)
% Defocus in diopters as a function of wavelength
%
% Syntax:
%   D = humanWaveDefocusS5(wave)
%
% Description:
%    This is a function fit to LCA data from the Chin et al. (2026) task.
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
% wavelength for human eye. Fit to data from subject S5 using
% ARCnlzLCA with random seed set at 1.

q1 = 1.2108;
q2 = 0.2612;    
q3 = 0.3601;

% This is the human defocus as a function of wavelength. This formula
% converts the wave in nanometers to wave in microns. D is in diopters.
Dref = q1 - (q2 ./ (waveRef * 1e-3 - q3));
D = q1 - (q2 ./ (wave * 1e-3 - q3))-Dref;

return
