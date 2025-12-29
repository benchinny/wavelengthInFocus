function negLL = psyfitWeibullNegLL(param,X,RC,aFix,bFix)

% function negLL = psyfitWeibullNegLL(param,X,RC,aFix,bFix)
%
%   example call:
%
% computes negative log likelihood
%
% param:          parameters for psyfitWeibullfunc
% X:              stimulus values
% RC:             responses
%                 1 -> correct
%                 0 -> incorrect
% aFix:           fixed sigma parameter value
% bFix:           fixed beta  parameter value
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% negLL:          negative log-likelihood


if  exist('aFix','var')    && ~isempty(aFix)    param(1) = aFix; end
if  exist('bFix','var')    && ~isempty(bFix)    param(2) = bFix; end

% COMPUTES LOG PROBABILITY OF EVERY TRIAL AND SUMS ACROSS TRIALS
negLL = -(sum(log(     psyfitWeibullfunc(X(RC==1),param(1),param(2),[]) )) + ...
          sum(log( 1 - psyfitWeibullfunc(X(RC==0),param(1),param(2),[]) )) );