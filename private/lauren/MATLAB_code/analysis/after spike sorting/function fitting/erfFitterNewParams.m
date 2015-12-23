function [estimatedParams projectionComplete error] = gaussCdfFitter(data, varargin)

% description: estimates parameters that give best fit of data to a cumulative gaussian function
%
% returns:
%   estimatedParams: [t w]
%   projectionsComplete: array with projectionsComplete(1,:) = set of stimulus amplitude values
%      ranging from first to last value given in data, and projectionsComplete(2,:) = expected response
%      probabilities based on curve fit
%   error: sum-squared residual values
%
% arguments:
%
% data 2-d array, with data(1,:) vector of stimulus amplitudes and data(2,:) vector of success rates
%
% parameter-value pairs:
% 
% 'tStart': define a starting value for the parameter t
% 'wStart': define a starting value for the parameter w
% 'makePlot': if set to 1, true, or 'true', generates plot of response curve
% 'lockedAmps': if not empty, data points in plot are black circles if analysis for that movie has
%   been "locked" (signified by a 1 value in "lockedAmps"), and grey circles otherwise (signified by
%   a 0 value in "lockedAmps")
%
% function has 2 parameters: w (slope factor) an t (threshold)
%   (response probability) = Phi[w(x/t - 1)] = 0.5*(1 + erf([w(x/t - 1)]/sqrt(2)))
%   where Phi is the cdf of the standard normal distribution (mean = 0, standard devation = 1), and
%   x represents the stimulus amplitude

% the actual curve form is identical to the previously-used erf, but the parameters have been
% repackaged:
%
%   (response rate) = 0.5*(1 + erf((stimulus amplitude)*a + b))
%   w = -sqrt(2)*b and t = -b/a
%
%
% Make a guess for initial estimate of w and t (start).  fminsearch
% minimizes the error returned from fitToErf by adjusting nlParams.  It returns the
% final value of nlParams.


%%

p = inputParser;

p.addRequired('data', @isnumeric)

p.addParamValue('tStart', 0.5, @isnumeric)
p.addParamValue('wStart', 5, @isnumeric)
p.addParamValue('makePlot', 0, @(x)any([x==0, x==1, strcmpi('true', x), strcmpi('false', x), islogical(x)]))
p.addParamValue('lockedAmps', [], @isnumeric)
p.addParamValue('setParam', 0, @isnumeric)

p.parse(data, aGuess, bGuess, varargin{:})


makePlot = p.Results.makePlot;
lockedAmps = p.Results.lockedAmps;
setParam = p.Results.setParam;
wStart = p.Results.wStart;
tStart = p.Results.tStart;

if strcmpi(makePlot, 'true')
    makePlot = 1;
elseif strcmpi(makePlot, 'false')
    makePlot = 0;
end

if setParam == 1 %hold parameter 1 (t) constant
    param1 = tStart;
    start = wStart;
elseif setParam == 2 %hold parameter 2 (w) constant
    param2 = wStart;
    start = tStart;
else
    start = [tStart wStart];
end


nAmps = size(data, 2);

% an anonymous function is used to specify which of the fitToErf arguments to search
% outputFcn = @(x,optimvalues,state) fitoutputfun(x,optimvalues,state,t,y,h);
% options = optimset('OutputFcn',outputFcn,'TolX',0.1);

if setParam == 1 %hold parameter 1 constant
    estimatedParam2 = fminsearch(@(param2)erfFun([param1 param2], data), start);
    estimatedParams = [param1 estimatedParam2];
elseif setParam == 2 %hold parameter 2 constant
    estimatedParam1 = fminsearch(@(param1)erfFun([param1 param2], data), start);
    estimatedParams = [estimatedParam1 param2];
else
    estimatedParams = fminsearch(@(nlParams)erfFun(nlParams, data), start);
end


a = estimatedParams(1);
b = estimatedParams(2);

xProj = data(1,1):0.001:data(1,nAmps);
projection = 0.5 + 0.5*erf(a*xProj+b);
projectionComplete = [xProj; projection]; 

error = norm(0.5 + 0.5*erf(a*data(1,:)+b) - data(2,:));

if makePlot
    figure
    hold on
    if any(data(1,:) < 0) %working in log scale
        title(['error = ', num2str(error), 10, 'threshold = ', num2str(exp(-estimatedParams(2)/estimatedParams(1)))])
    else %probably in linear scale
        title(['error = ', num2str(error), 10, 'threshold = ', num2str(-estimatedParams(2)/estimatedParams(1))])
    end
    if ~isempty(lockedAmps)
        for i = 1:size(data, 2)
            if lockedAmps(i)
                plot(data(1,i), data(2,i),'k.', 'markerSize', 10);
            else
                plot(data(1,i), data(2,i),'.', 'markerSize', 10, 'markerFaceColor', [1 1 1], 'markerEdgeColor', [0.5 0.5 0.5]);
            end
        end
    else
        plot(data(1,:), data(2,:),'k.', 'markerSize', 10, 'markerFaceColor', [1 1 1]);
    end
    plot(xProj, projection,'m-');
    set(gca, 'ylim', [0 1])
    xlabel('current amplitude (µA)')
    ylabel('response rate')
    hold off
    
end