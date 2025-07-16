function wave = humanWaveDefocusInvertARC(waveRef,D,subjNum)
% Defocus in diopters as a function of wavelength
%
% Syntax:
%   D = humanWaveDefocusInvert(wave)
%
% Description:
%    This is a function fit to the data from Bedford and Wyszecki and Wald
%    on human chromatic aberration.
%
%    This function contains examples of usage inline. To access, type 'edit
%    humanWaveDefocusInvert.m' into the Command Window.
%
% Inputs:
%    wave - Vector. Wavelength vector, in nanometers.
%
% Outputs:
%    D    - Vector. Defocus vector.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/11       Copyright ImagEval Consultants, LLC, 2011.
%    06/28/18  jnm  Formatting

% Example:
%{
    wave = 400:10:700;
    D = humanWaveDefocusInvert(wave);
    vcNewGraphWin;
    plot(wave, D);
    xlabel('Wave (nm)')
    ylabel('Diopters');
    grid on
%}

% Constants for formula to compute defocus in diopters (D) as a function of
% wavelength for human eye. From 
% Larry N. Thibos, Ming Ye, Xiaoxiao Zhang, and Arthur Bradley, 
% "The chromatic eye: a new reduced-eye model of ocular chromatic aberration 
% in humans," Appl. Opt. 31, 3594-3600 (1992)
if subjNum==10
    q1 = 1.9851;
    q2 = 0.4644;
    q3 = 0.3135;
elseif subjNum==1
    q1 = 2.0000;
    q2 = 0.5196;    
    q3 = 0.2717;
elseif subjNum==3
    q1 = 2.0000;
    q2 = 0.4701;    
    q3 = 0.2778;    
elseif subjNum==5
    q1= 1.1639;
    q2 = 0.2444;
    q3 = 0.3650;
elseif subjNum==9
    q1 = 1.0427;
    q2 = 0.1579;
    q3 = 0.3700;    
elseif subjNum==16
    q1= 2.0000;
    q2 = 0.6764;
    q3 = 0.2273;
elseif subjNum==17
    q1= 2.0000;
    q2 = 0.4996;
    q3 = 0.2842;   
elseif subjNum==18
    q1= 2.0000;
    q2 = 0.8258;
    q3 = 0.1967;    
elseif subjNum==20
    q1= 1.6553;
    q2 = 0.9900;
    q3 = 0.0337;    
else
    q1 = 1.7312;
    q2 = 0.63346;
    q3 = 0.21410;
end

% This is the human defocus as a function of wavelength. This formula
% converts the wave in nanometers to wave in microns. D is in diopters.
Dref = q1 - (q2 ./ (waveRef * 1e-3 - q3));
q1 = q1-Dref;
wave = ((q2./(q1-D))+q3)./0.001;

% plot(wave, D);
% grid;
% xlabel('Wavelength (nm)');
% ylabel('relative defocus (diopters)');

return
