function [estimatedParams projectionComplete error] = weibullFitter(data, aGuess, bGuess, cGuess, varargin)

% description: estimates parameters that give best fit of data to a Weibull distribution CDF 
%
% arguments:
%
% data 2-d array, with data(1,:) vector of stimulus amplitudes and data(2,:) vector of success rates
% a, b, c are the Weibull function parameters k, lambda, theta (scalars)
%    good start values are a = 5, b = 1, c = 0 (for linear scale data)
%
% parameter-value pairs:
% 
% 'makePlot': if set to 1, true, or 'true', generates plot of response curve
% 'lockedAmps': if not empty, data points in plot are black circles if analysis for that movie has
%   been "locked" (signified by a 1 value in "lockedAmps"), and grey circles otherwise (signified by
%   a 0 value in "lockedAmps")
% 'binData': bins data vertically (into response rate bins)--meant to be used for data from many
%   aligned curves (only affects plotting)
% 'setParam': constrains one of the parameters to equal the initial value. Use 1 to constain a, 2 to
%   constrain b, and 3 to constrain c
%
% function has 2 parameters:
%   (response rate) = 1 - exp(-((stimulus amplitude - c)/b)^a)
%   
%   this is the cdf of the Weibull distrubtion with an extra parameter c to allow for shifts along
%   the x axis
%
% Make a guess for initial estimate of a, b, c (start).  fminsearch
% minimizes the error returned from fitToErf by adjusting nlParams.  It returns the
% final value of nlParams.


%%

p = inputParser;

p.addRequired('data', @isnumeric)
p.addRequired('aGuess', @isnumeric)
p.addRequired('bGuess', @isnumeric)
p.addRequired('cGuess', @isnumeric)


p.addParamValue('makePlot', 0, @(x)any([x==0, x==1, strcmpi('true', x), strcmpi('false', x), islogical(x)]))
p.addParamValue('lockedAmps', [], @isnumeric)
p.addParamValue('setParam', 0, @isnumeric)
p.addParamValue('binData', 0, @(x)any([x==0, x==1, strcmpi('true', x), strcmpi('false', x), islogical(x)]))

p.parse(data, aGuess, bGuess, cGuess, varargin{:})


makePlot = p.Results.makePlot;
lockedAmps = p.Results.lockedAmps;
setParam = p.Results.setParam;
binData = p.Results.binData;

if strcmpi(makePlot, 'true')
    makePlot = 1;
elseif strcmpi(makePlot, 'false')
    makePlot = 0;
end

if setParam == 1 %hold parameter 1 constant
    lockedParam = aGuess;
    start = [bGuess cGuess];
elseif setParam == 2 %hold parameter 2 constant
    lockedParam = bGuess;
    start = [aGuess cGuess];
elseif setParam == 3
    lockedParam = cGuess;
    start = [aGuess bGuess];
else
    start = [aGuess bGuess cGuess];
end


options = optimset('MaxFunEvals', 1000);

if setParam == 1 %hold parameter 1 constant
    estimatedParam23 = fminsearch(@(nlParams)weibullFun(nlParams, data, lockedParam, 1), start, options);
    estimatedParams(1) = lockedParam;
    estimatedParams(2:3) = estimatedParam23;
elseif setParam == 2 %hold parameter 2 constant
    estimatedParam13 = fminsearch(@(nlParams)weibullFun(nlParams, data, lockedParam, 2), start, options);
    estimatedParams(2) = lockedParam;
    estimatedParams([1 3]) = estimatedParam13;
elseif setParam == 3 %hold parameter 3 constant
    estimatedParam12 = fminsearch(@(nlParams)weibullFun(nlParams, data, lockedParam, 3), start, options);
    estimatedParams(3) = lockedParam;
    estimatedParams(1:2) = estimatedParam12;
else
    estimatedParams = fminsearch(@(nlParams)weibullFun(nlParams, data), start, options);
end


a = estimatedParams(1);
b = estimatedParams(2);
c = estimatedParams(3);


threshold = b*((-log(0.5))^(1/a))+c;

nAmps = size(data, 2);
xProj = data(1,1):0.001:data(1,nAmps);

projection = zeros(size(xProj));
for i = 1:length(xProj)
    if (xProj(i) - c)/b >= 0
        projection(i) = 1 - exp(-((xProj(i)-c)/b).^a);
    else
        projection(i) = 0;
    end
end

projectionComplete = [xProj; projection]; 

p = zeros(1, size(data, 2));
for i = 1:size(data,2)
    if (data(1,i) - c)/b > 0
        p(i) = 1 - exp(-((data(1,i)-c)/b)^a);
    else
        p(i) = 0;
    end
end
error = norm(p - data(2,:));

%% binning data and taking meaning at different success rates

if binData

    dataBins = cell(50,1);
    dataBinMeans = zeros(50,1);
    dataBinStds = zeros(50,1);

    for i = 1:50
        for j = 1:size(data,2)
            if data(2,j) >= (i-1)*0.02 && data(2,j) < i*0.02
                dataBins{i} = [dataBins{i} data(1,j)];
            elseif data(2,j) == 1 && i == 50;
                dataBins{50} = [dataBins{50} data(1,j)];
            end
        end
        dataBinMeans(i) = mean(dataBins{i});
        dataBinStds(i) = std(dataBins{i});
    end
end

%% plotting

if makePlot
    figure
    hold on
    title(['error = ', num2str(error), 10, 'threshold = ', num2str(threshold), 10,...
        'k = ', num2str(a), ', lambda = ', num2str(b), ', theta = ' num2str(c)])
    if ~isempty(lockedAmps)
        for i = 1:size(data, 2)
            if lockedAmps(i)
                plot(data(1,i), data(2,i),'k.', 'markerSize', 10);
            else
                plot(data(1,i), data(2,i),'.', 'markerSize', 10, 'markerFaceColor', [1 1 1], 'markerEdgeColor', [0.5 0.5 0.5]);
            end
        end
    else
        if ~binData
            plot(data(1,:), data(2,:),'k.', 'markerSize', 10, 'markerFaceColor', [1 1 1]);
        else
            for i = 1:50;
                plot([dataBinMeans(i) - dataBinStds(i), dataBinMeans(i) + dataBinStds(i)], [i*0.02-0.01, i*0.02-0.01], 'b-')
                plot(dataBinMeans(i), i*0.02-0.01, 'bo')
            end
        end
    end
    
    plot(xProj, projection,'m-');
    set(gca, 'ylim', [0 1])
    xlabel('current amplitude (µA)')
    ylabel('response rate')
    hold off
end