function frqCpu = smpFrq(smpPerUnit,numSmp,bNonZeroOnly)

% function frqCpu = smpFrq(smpPerUnit,numSmp,bNonZeroOnly)
%
%   example call: cpd = smpFrq(128,128,1)
% 
% sampled frequencies given a sampling rate and a patch
% NOTE: frequencies should always be sampled at zero
%
% smpPerUnit:   sampling rate
% numSmp:       number of samples
% bNonZeroOnly: 1 -> returns only the non-zero frequencies
%               0 -> returns all frequencies (default)
% 
%               see smpPos.m
% 
%%%%%%%%%%%
% frqCpu:       sampled frequencies in cycles per unit

if ~exist('bNonZeroOnly','var') || isempty(bNonZeroOnly)
   bNonZeroOnly=0; 
end

% IF NUM PIX IS EVEN
if mod(numSmp,2) == 0 
    minPos = -smpPerUnit/2;
    maxPos =  smpPerUnit/2 - smpPerUnit/max(numSmp);
% IF NUM PIX IS ODD
elseif mod(numSmp,2) == 1
    minPos = -smpPerUnit/2;
    maxPos =  smpPerUnit/2;
else
    error(['smpFrq: WARNING! num pix must be integer valued. numPix = ' num2str(numSmp)]);
end
frqCpu   = linspace(minPos,maxPos,max(numSmp));

if bNonZeroOnly
   frqCpu = frqCpu(frqCpu>0); 
end