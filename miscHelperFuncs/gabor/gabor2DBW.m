function [g,G] = gabor2DBW(X,Y,x0,y0,frqCpd,ortDeg,phsDeg,BWoct,BWort,bNormalize,bPLOT)         

% function [g,G] = gabor2DBW(X,Y,x0,y0,frqCpd,ortDeg,phsDeg,BWoct,BWort,bNormalize,bPLOT)         
%
%   example call: % EVEN PHASE GABOR 
%                   gabor2DBW(smpPos(128,128),[],0,0,6,0,0,1.5,40*pi./180,1,1);         
%
%                 % ODD  PHASE GABOR (POSITIVE LOBE ON LEFT)
%                   gabor2DBW(smpPos(128,128),[],0,0,6,0,90,1.5,40*pi./180,1,1);
%
%                 % ODD  PHASE GABOR (POSITIVE LOBE ON RIGHT)
%                   gabor2DBW(smpPos(128,128),[],0,0,6,0,-90,1.5,40*pi./180,1,1);
%
% implementation based on Geisler_GaborEquations.pdf in /VisionNotes
%
% X:          X position  in visual degrees [nxm] matrix
%             if [   ]-> 1 deg patch size, 128 smpPerDeg
%             if [  n]-> 1 deg patch size,   n smpPerDeg
%             if [1xn]-> patch size specified by samples
% Y:          Y position  in visual degrees [nxm] matrix
% x0:         X position  of gabor center                      [ scalar ]
% y0:         Y position  of gabor center                      [ scalar ]      
% frqCpd:     frequency   in cycles per deg                    [ scalar ]
% ortDeg:     orientation in deg                               [ scalar ]                         
% phsDeg:     phase       in deg                               [ scalar ]       
% BWoct:      frequency   bandwidth in octaves                 [ scalar ]
% BWort:      orientation bandwidth in radians                 [ scalar ]
% bNormalize: 1 -> normalize to vector magnitude of 1 
%             0 -> don't
% bPLOT:      1 -> plot
%             0 -> not
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g:          gabor in space
% G:          gabor in frequency domain 
%             where G = fftshift(fft2(fftshift(g)))./sqrt(numel(g));

if ~exist('bPLOT','var') || isempty(bPLOT) bPLOT = 0; end
if isempty(X) || isempty(Y)
   if isempty(X)
   [X Y] = meshgrid(smpPos(128,128));    
   elseif isscalar(X)
   [X Y] = meshgrid(smpPos(X,X));
   elseif isvector(X)
   [X Y] = meshgrid(X);    
   elseif min(size(X)) > 1
   Y = X';    
   end
end
if isvector(X) && isvector(Y) 
   [X Y]=meshgrid(X,Y); 
end

% ORIENTATION IN RADIANS
ortRad = ortDeg.*pi./180;

% ROTATED POSITIONS: apply rotation matrix
Xp        =  (X-x0).*cos(ortRad) + (Y-y0).*sin(ortRad);
Yp        = -(X-x0).*sin(ortRad) + (Y-y0).*cos(ortRad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GAUSSIAN ENVELOPE STANDARD DEVIATIONS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIGMA IN SPATIAL FREQUENCY... BANDPASS DIRECTION (ORTHOGONAL TO GRATING)
sigmaXdeg = bandwidthOct2sigma(frqCpd, BWoct);
% SIGMA IN ORIENTATION...       LOWPASS  DIRECTION (PARALLEL   TO GRATING)  
sigmaYdeg = bandwidthOrt2sigma(frqCpd, BWort);

% GABOR
g = exp(-0.5.*(Xp./sigmaXdeg).^2)  .* ...
    exp(-0.5.*(Yp./sigmaYdeg).^2)  .* ...
    cos( (2.*pi.*frqCpd.*Xp) + phsDeg.*pi./180);
    
% NORMALIZE TO VECTOR MAGNITUDE OF 1
if bNormalize
    g = g./sqrt(sum(g(:).^2));
end

if nargout > 1 
    G = fftshift(fft2(ifftshift(g)))./sqrt(numel(g));
end
% COMPUTE BANDWIDTHS IN DEGREES
%
% ... DO THIS AT SOME POINT ...
%
%%
if bPLOT
   %%%%%%%%%%%%%%%%%%%%%%%
   % PLOT GABOR IN SPACE %
   %%%%%%%%%%%%%%%%%%%%%%%
   figure('position',[411  523  1014  531]);
   s1=subplot(1,2,1); hold on
   imagesc(X(1,:),Y(:,1)',g);
   % surf(X,Y,g,'edgecolor','none','facealpha',.8);
   formatFigure('X (deg)','Y (deg)','Space',0,0,18,14);
   axis square;
   grid on
   % view([-39    16]);
   view([0 90]);
%    rotate3d
    caxis(max(abs(minmax(g)))*[-1 1])
   axis xy
   xlim(minmax(X));
   try
   ylim(minmax(Y));
   end
   
   % COMPUTE FOURIER TRANSFORM
   G = fftshift(fft2(fftshift(g)))./sqrt(numel(g));
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   % PLOT GABOR IN FREQUENCY %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   subplot(1,2,2);
   
   % FIND NATIVE SAMPLING FREQUENCY IN EITHER DIRECTION
    smpPosPerUnitX = smpPos2smpPerUnit(X(1,:));
    smpPosPerUnitY = smpPos2smpPerUnit(Y(:,1));
   % NUMBER OF SAMPLES IN EITHER DIRECTION
    nx = size(X,2);
    ny = size(Y,1);
    % FREQUENCIES
    fxcpd = smpFrq(smpPosPerUnitX,nx);
    fycpd = smpFrq(smpPosPerUnitY,ny);
    %[U,V]=meshgrid(fxcpd,fycpd); 
    
   imagesc(fxcpd,fycpd,abs(G));
   formatFigure('U (cpd)','V (cpd)','Frequency',0,0,18,14);
   axis square;
   axis tight;
   axis xy
   view([0 90]);
   cb = colorbar;
   
   figure(gcf); suptitle([num2str(frqCpd,'%.2f') 'cpd' ...
                ', BW_{oct}=' num2str(BWoct,2) ...
                ', BW_{\theta}=' num2str(BWort.*180./pi,2)   ...
                ', \Theta=' num2str(ortDeg,3) ...
                ', \Phi=' num2str(phsDeg,3)],22);
   set(cb,'position',[ 0.9159    0.1714    0.0247    0.6405]);
   subplot(s1);
   colormap(cmapBWR(256));
end
