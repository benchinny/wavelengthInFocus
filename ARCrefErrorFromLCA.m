function deltaR = ARCrefErrorFromLCA(lambda)

% function deltaR = ARCrefErrorFromLCA(lambda)
%
% example call: deltaR = ARCrefErrorFromLCA(589)
%
% for computing defocus as a function of wavelength. Implements equation
% from Thibos et al. (1992) and Roorda et al. (2023)

n_ref = 1.33413;
r = 0.00555;
n_D = 1.333;
a = 1.31848;
b = 0.006662;
c = 0.1292;

deltaR = (n_ref-a-b./(lambda-c))./(r*n_D);

end