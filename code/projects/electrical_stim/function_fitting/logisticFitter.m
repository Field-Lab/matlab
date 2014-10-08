function [estimatedParams projectionComplete error] = logisticFitter(data, aGuess, bGuess, varargin)

% description: estimates parameters that give best fit of data to an error function with a lower
% asymptote at 0 and an upper asymptote at 1
%
% arguments:
%
% data 2-d array, with data(1,:) vector of stimulus amplitudes and data(2,:) vector of success rates
% a and b are the error function parameters (scalars)
%
% parameter-value pairs:
% 
% 'makePlot': if set to 1, true, or 'true', generates plot of response curve
% 'lockedAmps': if not empty, data points in plot are black circles if analysis for that movie has
%   been "locked" (signified by a 1 value in "lockedAmps"), and grey circles otherwise (signified by
%   a 0 value in "lockedAmps")
%
% function has 2 parameters:
%   (response rate) = 0.5 + 0.5*erf((stimulus amplitude)*a + b)
%   
%   this is an error function constrained to run from 0 (-infinity end) to 1 (+infinity end)
%
% Make a guess for initial estimate of a and b (start).  fminsearch
% minimizes the error returned from fitToErf by adjusting nlParams.  It returns the
% final value of nlParams.


%%

p = inputParser;

p.addRequired('data', @isnumeric)
p.addRequired('aGuess', @isnumeric)
p.addRequired('bGuess', @isnumeric)


p.addParamValue('makePlot', 0, @(x)any([x==0, x==1, strcmpi('true', x), strcmpi('false', x), islogical(x)]))
p.addParamValue('lockedAmps', [], @isnumeric)
p.addParamValue('setParam', 0, @isnumeric)

p.parse(data, aGuess, bGuess, varargin{:})


makePlot = p.Results.makePlot;
lockedAmps = p.Results.lockedAmps;
setParam = p.Results.setParam;

if strcmpi(makePlot, 'true')
    makePlot = 1;
elseif strcmpi(makePlot, 'false')
    makePlot = 0;
end

if setParam == 1 %hold parameter 1 constant
    param1 = aGuess;
    start = bGuess;
elseif setParam == 2 %hold parameter 2 constant
    param2 = bGuess;
    start = aGuess;
else
    start = [aGuess bGuess];
end

nAmps = size(data, 2);

% an anonymous function is used to specify which of the fitToErf arguments to search
% outputFcn = @(x,optimvalues,state) fitoutputfun(x,optimvalues,state,t,y,h);
% options = optimset('OutputFcn',outputFcn,'TolX',0.1);

if setParam == 1 %hold parameter 1 constant
    estimatedParam2 = fminsearch(@(param2)logisticFun([param1 param2], data), start);
    estimatedParams = [param1 estimatedParam2];
elseif setParam == 2 %hold parameter 2 constant
    estimatedParam1 = fminsearch(@(param1)logisticFun([param1 param2], data), start);
    estimatedParams = [estimatedParam1 param2];
else
    estimatedParams = fminsearch(@(nlParams)logisticFun(nlParams, data), start);
end


a = estimatedParams(1);
b = estimatedParams(2);

xProj = data(1,1):0.001:data(1,nAmps);
%projection = 1./(1 + exp(-(((xProj-b)./a).^2)));
projectedData = zeros(1, size(data, 2));
projection = zeros(1, length(xProj));
for i = 1:length(xProj)
    if (xProj(i)-b)/a>0
        projection(i) = 1./(1 + exp(-(((xProj(i)-b)./a).^2)));
    else
        projection(i) = 1./(1 + exp((((xProj(i)-b)./a).^2)));
    end
end

for i = 1:length(projectedData)
    if (data(1,i)-b)/a>0
        projectedData(i) = 1./(1 + exp(-(((data(1,i)-b)./a).^2)));
    else
        projectedData(i) = 1./(1 + exp((((data(1,i)-b)./a).^2)));
    end
end
    
    
projectionComplete = [xProj; projection]; 

%error = norm(1./(1 + exp(-(((data(1,:)-b)./a).^2))) - data(2,:));
%error = norm(1./(1 + exp(-(((data(1,:).^2-b)./a).^2))) - data(2,:));
error = norm(projectedData - squeeze(data(2,:)));


if makePlot
    figure
    hold on
    title(['error = ', num2str(error), 10, 'threshold = ', num2str(-estimatedParams(2)/estimatedParams(1))])
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