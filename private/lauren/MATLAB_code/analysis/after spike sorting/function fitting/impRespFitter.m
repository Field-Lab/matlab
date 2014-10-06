function [estParams alpha_const FWHM ttp error] = impRespFitter(data, t0Guess, tauGuess, n, varargin)
%function [estParams alpha error] = impRespFitter(data, tauGuess, n, varargin)


% description: estimates parameters that give best fit of data to the
% impulse response of a set of low-pass filters

% arguments:
%
% data: an array with dimenstions 2 x N, with data(1,:) = x-values and
% data(2,:) = y-values
%  *****x-values must be linearly increasing!
%
%
% t0 reflects latency, tau reflects width, alpha reflects amplitude, and n
% is the number of filters in the sequence
%
% alpha is constrained to give approximately the same area under the curve
% as under the histogram being fit
%
% FWHM = full width at half maximum
% ttp = time to peak
%
% parameter-value pairs:
% 
% 'makePlot': if set to 1, true, or 'true', generates plot of response
% curve
%
%
% Make a guess for initial estimate of t0 and tau (start).  fminsearch
% minimizes the error returned from impRespFun by adjusting nlParams.  It returns the
% final value of nlParams.

%%

p = inputParser;

p.addRequired('data', @isnumeric)
p.addRequired('t0Guess', @isnumeric)
p.addRequired('tauGuess', @isnumeric)
p.addRequired('n', @isnumeric)



p.addParamValue('makePlot', 0, @(x)any([x==0, x==1, islogical(x)]))
p.addParamValue('plotAxes', [], @ishandle)

p.parse(data, t0Guess, tauGuess, n, varargin{:})
%p.parse(data, tauGuess, n, varargin{:})

makePlot = p.Results.makePlot;
plotAxes = p.Results.plotAxes;

start = [t0Guess tauGuess];
%start = tauGuess

% an anonymous function is used to specify which of the fitToErf arguments to search
% outputFcn = @(x,optimvalues,state) fitoutputfun(x,optimvalues,state,t,y,h);
options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, 'FunValCheck', 'on');

estParams = fminsearch(@(nlParams)impRespFun(nlParams, n, data), start, options);
estParamsConst = fminsearch(@(nlParams)impRespFunConstArea(nlParams, n, data), start, options);


t0 = estParams(1);
estParams(2) = abs(estParams(2));
tau = estParams(2);

t0_const = estParamsConst(1);
estParamsConst(2) = abs(estParamsConst(2));
tau_const = estParamsConst(2);

%determine alpha for unconstrained area case
t_sc = (data(1,:) - t0)/tau; %scaled/shifted time values
p = exp(-n*(t_sc-1)).*t_sc.^n;
p(t_sc<0) = 0; %zero out values corresonding to t_xc < 0

% %find alpha that gives lowest error
A = p';
b = data(2,:)';
alpha = A\b;

%find alpha that gives equal area under curve and histogram *** assumes
%that x-values in data are linearly increasing!

%subsample curve fit (divide each interval into 100 intervals)
xVals = linspace(data(1,1), data(1,end), 100*(size(data,2)-1)+1);
xVals = (xVals - t0_const)/tau_const; %scale according to nonlinear params
pSubSamp = exp(-n*(xVals-1)).*xVals.^n;
pSubSamp(xVals<0) = 0;
pSubSampArea = sum(pSubSamp)/100; %scaled by 100 to account for subsampling

alpha_const = sum(data(2,:))/pSubSampArea;

%calculate squared error
error = norm(alpha*p - data(2,:));



%numerically calculate FWHF (for unconstrained fit)
tProj = 0:0.0001:data(1,end); %determines precision of FWHM estimate

%unconstrained fit
t_sc = (tProj - t0)/tau; %scaled/shifted time values
p = exp(-n*(t_sc-1)).*t_sc.^n;
p(t_sc<0) = 0; %zero out values corresonding to t_sc < 0

HFvals = [find(p>=0.5, 1, 'first') find(p>=0.5, 1, 'last')];

FWHM = tProj(HFvals(2))-tProj(HFvals(1));
ttp = tProj(find(p == max(p), 1));

if makePlot
%     tProj = 0:0.001:data(1,end);
%     
%     %unconstrained fit
%     t_sc = (tProj - t0)/tau; %scaled/shifted time values    
%     p = exp(-n*(t_sc-1)).*t_sc.^n;
%     p(t_sc<0) = 0; %zero out values corresonding to t_xc < 0
    
    %constrained fit
    t_sc_const = (tProj-t0_const)/tau_const;
    p_const = exp(-n*(t_sc_const-1)).*t_sc_const.^n;
    p_const(t_sc_const<0) = 0;
    
    if ~isempty(plotAxes)
        axes(plotAxes); hold on
    else
        figure('position', [400 400 600 400])
        hold on
    end
    %plot(data(1,:), data(2,:),'k.', 'markerSize', 5, 'markerFaceColor', [1 1 1]);
    plot(tProj, alpha*p,'m-');
    plot(tProj(HFvals(1)) + [0 FWHM], 0.5*alpha*[1 1], 'm-')
    plot(tProj, alpha_const*p_const, 'b-')
    text(0.05, alpha, ['FWHM: ' num2str(FWHM)])
    text(0.05, alpha*0.9, ['time to peak: ' num2str(ttp)])
    
    hold off
end

