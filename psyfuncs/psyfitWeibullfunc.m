function [PC,T] = psyfitWeibullfunc(X,aFit,bFit,PCcrt)

% function [PC,T] = psyfitWeibullfunc(X,aFit,bFit,PCcrt)
%
% plot Weibull fit to psychometric data
%
% X:         stimulus values                             [ 1 x nTrials ]
% aFit:      alpha parameter in Weibull function         [1  x     1]
% bFit:      beta of Weibull function
% PCcrt:     criterion performance corresponding to threshold
%%%%%%%%%%%%%%%%%%%%%%
% PC:        percent correct                              [ 
% T:         threshold corresponding to criterion performance

if ~exist('PCcrt','var') || isempty(PCcrt), PCcrt = 0.85; end

% TYPICAL EQUATION FOR A WEIBULL FUNCTION WITH A LAPSE RATE OF 0.5, E.G.
% SEE MORTENSEN (2002) IN VISION RESEARCH.
PC0 = 0.5;
PC = 1-(1-PC0).*exp(-(X./aFit).^bFit);

% SIGNAL AT THRESHOLD (STANDARD EQUATION)
T  = abs(-aFit.*log((1-PCcrt)./(1-PC0)).^(1./bFit));
