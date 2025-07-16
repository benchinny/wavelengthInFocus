function [gRGB,G] = ARC2Dgabor(x,y,x0,y0,frqCpd,A,ortDeg,phsDeg,sigmaX,sigmaY,rgb,bFixSigmaX,bFixSigmaY,bNormalize,bPLOT)

% function [gRGB,G] = ARC2Dgabor(x,y,x0,y0,frqCpd,A,ortDeg,phsDeg,sigmaX,sigmaY,rgb,bFixSigmaX,bFixSigmaY,bNormalize,bPLOT)    
%
%   example call: [gRGB G] = ARC2Dgabor(smpPos(128,128),[],0,0,[9],[1],[-15],[0],0.187,0.187,[0.56 0 1.00],1,1,0,1);
%
% implementation based on Geisler_GaborEquations.pdf in /VisionNotes
%
% x:          x position  in visual degrees matrix
%             []    -> 1 deg patch size, 128 smpPerDeg
%             [n]   -> 1 deg patch size,   n smpPerDeg
%             [1xn] -> x(end)-x(1) patch size, length(x) samples
% y:          y position  in visual degrees matrix       
% x0:         x position  of gabor center                      [  scalar  ] 
% y0:         y position  of gabor center                      [  scalar  ]
% frqCpd:     frequency   in cycles per deg                    [1 x nComp ]
% A     :     amplitude 
%             [scalar] -> assigns same amplitude to all components
%             [1 x nComp] -> unique amplitude for each component
% ortDeg:     orientation in deg                               
%             [scalar] -> assigns same orientation to all components
%             [1 x nComp] -> unique orientation for each component
% phsDeg:     phase       in deg                                
%             [scalar] -> assigns same phase to all components
%             [1 x nComp] -> unique phase for each component
% sigmaX:     horizontal width of Gaussian envelope         
%             [scalar] -> assigns same bandwidth to all components
%             [1 x nComp] -> unique bandwidth for each component
% sigmaY:     vertical width of Gaussian envelope
%             [scalar] -> assigns same bandwidth to all components
%             [1 x nComp] -> unique bandwidth for each component
% rgb   :     color in rgb values
%             [1 x 3]
% bFixSigmaX: 1 -> fix sigmaX according to fundamental freq oct bandwidth
%             0 -> don't
% bFixSigmaY: 1 -> fix sigmaY according to fundamental freq ornt bandwidth
%             0 -> don't
% bNormalize: 1 -> normalize to vector magnitude of 1
%             0 -> don't
% bPLOT:      1 -> plot
%             0 -> not
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g:          gabor
% G:          gabor in frequency domain 
%             where G = fftshift(fft2(fftshift(g)))./sqrt(numel(g));

if ~exist('bPLOT','var') || isempty(bPLOT) bPLOT = 0; end
if isempty(x) || isempty(y)
   if isempty(x)
   [x, y] = meshgrid(smpPos(128,128));    
   elseif isscalar(x)
   [x, y] = meshgrid(smpPos(x,x));
   elseif isvector(x)
   [x, y] = meshgrid(x);    
   elseif min(size(x)) > 1
   y = x';    
   end
end
if isvector(x) && isvector(y) 
   [x, y]=meshgrid(x,y); 
end
% DEFAULT PARAMETER VALUES
if isempty(x0); x0 = 0; end
if isempty(y0); y0 = 0; end
if isempty(frqCpd); error('ARC2Dgabor: SPECIFY frqCpd!'); end
if isempty(A); A = ones([1 length(frqCpd)]); end
if isempty(ortDeg); ortDeg = zeros([1 length(frqCpd)]); end
if isempty(phsDeg); phsDeg = zeros([1 length(frqCpd)]); end

% IF ortDeg IS SCALAR
if length(ortDeg)==1
   ortDeg = ortDeg.*ones([1 length(frqCpd)]); 
end

% IF phsDeg IS SCALAR
if length(phsDeg)==1
   phsDeg = phsDeg.*ones([1 length(frqCpd)]); 
end

% IF A IS SCALAR
if length(A)==1
   A = A.*ones([1 length(frqCpd)]); 
end

% MAKE SURE PARAMETER VECTORS ARE SAME LENGTH
if ~(length(frqCpd)==length(ortDeg) && length(ortDeg)==length(phsDeg) && length(phsDeg)==length(A))
   error(['ARC2Dgabor: NUMBER OF ELEMENTS MUST BE SAME IN frqCpd, ortDeg,' ...
          ', phsDeg, and A!']); 
end

% THROWS ERROR BECAUSE IF FIXING sigmaX, CANNOT SPECIFY MORE THAN ONE
% sigmaX
if bFixSigmaX == 1 && length(sigmaX)>1
   error(['ARC2Dgabor: IF bFixSigmaX == 1, THEN sigmaX SHOULD CONTAIN' ...
          ' ONE ELEMENT ONLY--FOR THE FUNDAMENTAL FREQ!']); 
elseif bFixSigmaX~=1 && length(sigmaX)~=length(frqCpd)
   error(['ARC2Dgabor: sigmaX SHOULD CONTAIN THE SAME NUMBER OF ELEMENTS' ...
          ' AS THE OTHER STIMULUS PARAMETERS (E.G. FRQCPD) IF bFixSigmaX~=1']);
end

% % THROWS ERROR BECAUSE IF FIXING sigmaX, CANNOT SPECIFY MORE THAN ONE
% sigmaY
if bFixSigmaY == 1 && length(sigmaY)>1
      error(['ARC2Dgabor: IF bFixSigmaY == 1, THEN sigmaY SHOULD CONTAIN' ...
             ' ONE ELEMENT ONLY--FOR THE FUNDAMENTAL FREQ!']); 
elseif bFixSigmaY~=1 && length(sigmaY)~=length(frqCpd)
      error(['ARC2Dgabor: sigmaY SHOULD CONTAIN THE SAME NUMBER OF ELEMENTS' ...
             ' AS THE OTHER STIMULUS PARAMETERS (E.G. FRQCPD) IF bFixSigmaY~=1']);
end

if bFixSigmaX == 1 % IF FIXING sigmaX
   % CONVERT COMMON sigmaX VALUE TO OCTAVE BANDWIDTH
   BWoct = sigma2bandwidthOct(frqCpd,sigmaX); 
end

if bFixSigmaY == 1 % IF FIXING sigmaY
   % CONVERT COMMON sigmaY VALUE TO ORIENTATION BANDWIDTH
   BWort = sigma2bandwidthOrt(frqCpd,sigmaY); 
end

for i = 1:length(frqCpd) % FOR EACH FREQUENCY
    % CREATE GABOR
    gCmp(:,:,i) = A(i).*gabor2DBW(x,y,x0,y0,frqCpd(i),ortDeg(i),phsDeg(i),BWoct(i),BWort(i),0,0);
end

% ADD COMPONENTS
g = sum(gCmp,3);

% NORMALIZE TO MAGNITUDE OF 1
if bNormalize
    g = g./sqrt(sum(g(:).^2));
end

if nargout > 1 
    G = fftshift(fft2(ifftshift(g)))./sqrt(numel(g));
end

% MEAN LUMINANCE
DC = 0.5;
% TURN CONTRAST IMAGE INTO LUMINANCE IMAGE
gLum = g.*DC + DC;
% CROP CIRCLE?
bCircCrop = true;
if bCircCrop
    indCrop = sqrt(x.^2+y.^2) < 0.49;
    gLum(~indCrop) = 0;
end

% CONVERT INTO COLOR IMAGE
gRGB = [];
gRGB(:,:,1) = gLum.*rgb(1);
gRGB(:,:,2) = gLum.*rgb(2);
gRGB(:,:,3) = gLum.*rgb(3);

if bPLOT
   %%%%%%%%%%%%%%%%%%%%%%%
   % PLOT GABOR IN SPACE %
   %%%%%%%%%%%%%%%%%%%%%%%
   figure('position',[411  523  1014  531]);
   s1=subplot(1,2,1); hold on
   imagesc(x(1,:),y(:,1)',gRGB);
   % surf(X,Y,g,'edgecolor','none','facealpha',.8);
   formatFigure('X (deg)','Y (deg)','Space',0,0,18,14);
   axis square;
   grid on
   % view([-39    16]);
   view([0 90]);
%  rotate3d
   caxis(max(abs(minmax(g)))*[-1 1])
   axis xy
   xlim(minmax(x));
   ylim(minmax(y));
   
   % COMPUTE FOURIER TRANSFORM
   G = fftshift(fft2(fftshift(g)))./sqrt(numel(g));
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   % PLOT GABOR IN FREQUENCY %
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   subplot(1,2,2);
   
   % FIND NATIVE SAMPLING FREQUENCY IN EITHER DIRECTION
    smpPosPerUnitX = smpPos2smpPerUnit(x(1,:));
    smpPosPerUnitY = smpPos2smpPerUnit(y(:,1));
   % NUMBER OF SAMPLES IN EITHER DIRECTION
    nx = size(x,2);
    ny = size(y,1);
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
   
    figure(gcf); suptitle(['BW_{oct}=' '[' num2str(BWoct,2) ']' ...
                 ', BW_{\theta}=' '[' num2str(BWort.*180./pi,2) ']'   ...
                 ', \Theta=' '[' num2str(ortDeg,3) ']' ...
                 ', \Phi=' '[' num2str(phsDeg,3) ']'],22);
    set(cb,'position',[ 0.9159    0.1714    0.0247    0.6405]);
    subplot(s1);
%    colormap(cmapBWR(256));
end

end