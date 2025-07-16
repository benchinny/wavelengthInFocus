%%

q1 = 1.7312;
% q2 = 0.63346;
% q3 = 0.21410;
q2 = 0.01:0.01:0.9;
q3 = 0.37;
wave = 450:630;

% This is the human defocus as a function of wavelength. This formula
% converts the wave in nanometers to wave in microns. D is in diopters.
figure;
hold on;
for i = 1:length(q2)
    D = q1 - (q2(i) ./ (wave * 1e-3 - q3));
    plot(wave,D);
end

%%

clear;