function [aFit,bFit,Tfit,PCdta,PCfit,negLL] = psyfitWeibull(X,RC,aFix,bFix,PCcrt,bPLOT,xLbl,yLbl,color,shape,figh)

% function [aFit,bFit,Tfit,PCdta,PCfit,negLL] = psyfitWeibull(X,RC,aFix,bFix,PCcrt,bPLOT,xLbl,yLbl,color,shape,figh)
%
% fit generalized weibull function to data (modification of
% psyfitgengauss.m in Johannes Burge's lab toolbox by Ben Chin)
%
% X:         stimulus values                 [nTrl x 1]
% RC:        subject responses (0 or 1)      [nTrl x 1]
%            coded as correct vs incorrect
%            1 -> correct
%            0 -> incorrect
% aFix:      fixed value of alpha parameter
% bFix:      fixed value of beta           
% PCcrt:     criterion performance corresponding to threshold
% bPLOT:     1 -> plot
%            0 -> not
% xLbl:      plot x-label
% yLbl:      plot y-label
% color:     plot color
% shape:     plot shape
% figh:      plot figure handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% aFit:     alpha     fit
% bFit:     beta      fit
% Tfit:     threshold fit
% PCdta:    percent correct for the raw data
% PCfit:    percent correct for the fit at each X
% negLL:    negative log-likelihood of fit

% INPUT HANDLING
if ~exist('aFix','var')    || isempty(aFix);      aFix    =  []; end
if ~exist('bFix','var')    || isempty(bFix);      bFix    =  []; end
if ~exist('PCcrt','var')   || isempty(PCcrt),     PCcrt   =0.85; end
if ~exist('bPLOT','var')   || isempty(bPLOT),     bPLOT   =   0; end
if ~exist('xLbl','var');                          xLbl    =  []; end
if ~exist('yLbl','var');                          yLbl    =  []; end
if ~exist('color','var');                         color   = 'k'; end
if ~exist('shape','var');                         shape   = 'o'; end
if ~exist('figh','var');                          figh    =  []; end

% INPUT CHECKING
if size(X,2)  ~= 1, X  = X(:);    end
if size(RC,2) ~= 1, RC = RC(:); end
if length(X)  ~= length(RC), error(['psyfitWeibull: WARNING! X and RC sizes do not match. Fix it!']); end

% SET LOWER AND UPPER BOUNDS ON PARAMETERS
pLB     = [0.02.*(max(X(:))-min(X(:))) 0.35];
pUB     = [2.00.*(max(X(:))-min(X(:))) 3.00];

% SET FMINCON OPTIONS
minFuncType = 'fmincon';
if strcmp(minFuncType,'fmincon')
    opts             = optimset('fmincon');
    opts.Algorithm   = 'active-set';
    opts.LargeScale  = 'off';
    opts.UseParallel = 'never';
    opts.Display     = 'none';
    opts.MaxIter     = 500;
elseif strcmp(minFuncType,'fminsearch')
    opts             = optimset('fminsearch');
    opts.UseParallel = 'never';
    opts.Display     = 'off';
    opts.MaxIter     = 500;
end

% SET INITIAL PARAMETER VALUES
a0  = aFix;
b0  = bFix;
if isempty(a0); a0 = diff(minmaxLocal(abs(X)))./6;      a0 = a0  + .1.*a0.*randn; end
if isempty(b0); b0 = 1;                                 b0 = b0  + .1.*b0.*randn; end
p0 = [a0 b0];

% MINIMIZE NEGATIVE LOG-LIKELIHOOD
if strcmp(minFuncType,'fmincon')
    [pFit,negLL] = fmincon(   @(p) psyfitWeibullNegLL(p,X,RC,aFix,bFix),p0,[],[],[],[],pLB,pUB,[],opts);
elseif strcmp(minFuncType,'fminsearch')
    [pFit,negLL] = fminsearch(@(p) psyfitWeibullNegLL(p,X,RC,aFix,bFix),p0,opts);
end

% FINAL FIT PARAMETERS
if isempty(aFix); aFit = pFit(1); else; aFit = aFix; end
if isempty(bFix); bFit = pFit(2); else; bFit = bFix; end

% FIT THE FUNCTION
Xunq = unique(X);

% MAKE FIT WITH THRESHOLDS
[PCfit,Tfit] = psyfitWeibullfunc(Xunq,aFit,bFit,PCcrt);

% RAW DATA- COMPUTE PERCENT COMPARISON CHOSEN
[PCdta,~,Xunq] = psyPercentChosen(zeros(size(X)),X,RC);

%%%%%%%%%%%%%%%%
% PLOT RESULTS %
%%%%%%%%%%%%%%%%
if bPLOT
    %% OPEN FIGURE
    if ~exist('figh','var') || isempty(figh)
    figure('position',[680   634   384   406]); hold on
    else
	figure(figh);  hold on
    end
    % PLOT FIT (IN HI-RES)
    XcmpPlt = linspace(1.5.*min(X-mean(X))+mean(X),1.5.*max(X-mean(X))+mean(X),201);
    [PCplt,T]=psyfitWeibullfunc(XcmpPlt,aFit,bFit,PCcrt); hold on;
    plot(XcmpPlt,PCplt,'color',color,'linewidth',1.5);
    % PLOT DATA
    if strcmp(shape,'-') || strcmp(shape,'--') shape = 'o'; end
    plot(Xunq,PCdta,shape,'color',color,'linewidth',2,'markersize',15,'markerface','w');
    % WRITE STUFF TO SCREEN
    if isempty(figh)
    writeText(1-.1,.1,{['n=' num2str(numel(RC))]},'ratio',18,'right')
    end
    formatFigure([xLbl],[yLbl],['T=' num2str(T,'%.2f') ': ,\sigma=' num2str(aFit,'%1.2f') ',\beta=' num2str(bFit,'%1.2f')]);
    xlim(minmaxLocal(X)+[-.1 .1]); ylim([0 1])
    axis square
end
