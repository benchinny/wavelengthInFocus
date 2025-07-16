function wave = humanWaveDefocusInvertSteep(D)

% D = 2*((humanWaveDefocus(wave)+1.7).^0.3)-2.2;

wave = humanWaveDefocusInvert(((D+8)./8).^(1/0.05)-1.7);

% lambda = -log(1./D-1);
% wave = 14.*lambda+550;

end