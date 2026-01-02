function err = humanWaveDefocusObjFunc(wave,D,p)

% objective function for fitting LCA curves. Takes in wavelength, defocus,
% and parameters from Thibos (1992) model. Returns the error between the
% model's defocus values and the input defocus values.

q1 = p(1);
q2 = p(2);
q3 = p(3);

% This is the human defocus as a function of wavelength. This formula
% converts the wave in nanometers to wave in microns. D is in diopters.
Dpred = q1 - (q2 ./ (wave * 1e-3 - q3));
err = sqrt(mean((D-Dpred).^2));

return
