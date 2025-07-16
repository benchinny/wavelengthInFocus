function D = humanWaveDefocusSteep(wave)

D = 8*((humanWaveDefocus(wave)+1.7).^0.05)-8;

% D = 1./(1+exp(-(wave-550)./14));

end