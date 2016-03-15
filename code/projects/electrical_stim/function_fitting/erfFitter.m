function [estimatedParams projectionComplete error] = erfFitter(data, aGuess, bGuess, varargin)

% description: estimates parameters that give best fit of data to an error function with a lower
% asymptote at 0 and an upper asymptote at 1 by maximizing the log
% likelihood function
%
% arguments:
%
% data: an array with dimenstions 3 x nAmps, with data(1,:) vector of stimulus amplitudes, data(2,:) vector of success rates, and
%       data(3,:) vector of number of trials
%
%
%  *** maximum likelihood fitting is used to account for differences in certainty in
%  probability estimates at different amplitudes (based on Watson 1979 in
%  Vision Research 19:515-522)
% 
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
%   (response rate) = 0.5*(1 + erf((stimulus amplitude)*a + b))
%   this is an error function constrained to run from 0 (-infinity end) to 1 (+infinity end)
%
%   note: cumulative gaussian: Phi(y) = 0.5*(1 + erf(y/sqrt(2)))
%     if y is expressed as w(x/t - 1), then w = -sqrt(2)*b and t = -b/a
%     standard deviation (sigma) of cumulative gaussian = 1/(a*sqrt(2)) = t/w
%
%
% Make a guess for initial estimate of a and b (start).  fminsearch
% minimizes the error returned from fitToErf by adjusting nlParams.  It returns the
% final value of nlParams.
%
% Oct. 2011: changed from minimizing squared residuals to maximizing likelihood
% Dec. 2011: added 'useRobustFitting' option


%%

p = inputParser;

p.addRequired('data', @isnumeric)
p.addRequired('aGuess', @isnumeric)
p.addRequired('bGuess', @isnumeric)


p.addParamValue('makePlot', 0, @(x)any([x==0, x==1, strcmpi('true', x), strcmpi('false', x), islogical(x)]))
p.addParamValue('lockedAmps', [], @isnumeric)
p.addParamValue('setParam', 0, @isnumeric)
p.addParamValue('useRobustFitting', true, @islogical) %limit data points used for fitting to those that lie between 10% and 90%,
%(based on curve fit), plus one data point above and below this range

p.parse(data, aGuess, bGuess, varargin{:})


makePlot = p.Results.makePlot;
lockedAmps = p.Results.lockedAmps;
setParam = p.Results.setParam;
useRobustFitting = p.Results.useRobustFitting;

%limits used by robust fitting
lowerProbLim = 0.1;
upperProbLim = 0.9;


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
options = optimset('MaxFunEvals', 5000, 'MaxIter', 5000, 'FunValCheck', 'on');

if size(data,1) == 3 %maximize log-likelihood
    if setParam == 1 %hold parameter 1 constant
        estimatedParam2 = fminsearch(@(param2)erfLikelihood([param1 param2], data), start, options);
        estimatedParams = [param1 estimatedParam2];
    elseif setParam == 2 %hold parameter 2 constant
        estimatedParam1 = fminsearch(@(param1)erfLikelihood([param1 param2], data), start, options);
        estimatedParams = [estimatedParam1 param2];
    else
        estimatedParams = fminsearch(@(nlParams)erfLikelihood(nlParams, data), start, options);
    end
else
    errordlg('Unexpected dimensionality of argument "data": Make sure you have provided the number of trials for each amplitude!')
end


if useRobustFitting 
    %limit data points used for fitting to those that lie between 0.1 and 0.9 response probability 
    % (based on curve fit), plus one data point above and below this range
    exitRobustFitting = false;
        
    % determine which data points lie between lowerProbLim and
    % upperProb Lim (based on curve fit), plus one data point above and
    % below this range
    
    projectedProbs = 0.5 + 0.5*erf(estimatedParams(1)*data(1,:)+estimatedParams(2));
    
    if ~any(projectedProbs > lowerProbLim) || ~any(projectedProbs < upperProbLim)
        exitRobustFitting = true;
    end
            
    aboveMin = projectedProbs > lowerProbLim;
    belowMax = projectedProbs < upperProbLim;
    
    minMovInRangeInd = find(aboveMin,1);
    maxMovInRangeInd = find(belowMax,1,'last');
    
    firstMovInd = max([minMovInRangeInd-1, 1]);
    lastMovInd = min([maxMovInRangeInd+1, length(projectedProbs)]);

    iter = 0;
    repeatedRangeInd = 0;
    while ~exitRobustFitting
        iter = iter + 1;
        %disp(['robust fitter iteration ' num2str(iter)])
        
        %checks if current limits have been used in the past
        for ii = 1:iter-1
            if firstMovInd == firstMovHist(ii) && lastMovInd == lastMovHist(ii)
                repeatedRangeInd = ii;
                break
            end
        end
        
        if repeatedRangeInd %looping over same sets of limits
            %choose range that maximizes number of included data points
            firstMovInd = min(firstMovHist(repeatedRangeInd:end));
            lastMovInd = max(lastMovHist(repeatedRangeInd:end));
            exitRobustFitting = true; %refit one final time using this range
        end
        
        firstMovHist(iter) = firstMovInd; %#ok<AGROW>
        lastMovHist(iter) = lastMovInd; %#ok<AGROW>
        
        dataInRange = data(:,firstMovInd:lastMovInd);

        %redo fitting using only data within range
        if size(data,1) == 3 %maximize log-likelihood
            if setParam == 1 %hold parameter 1 constant
                estimatedParam2 = fminsearch(@(param2)erfLikelihood([param1 param2], dataInRange), estimatedParams(2), options);
                estimatedParams = [param1 estimatedParam2];
            elseif setParam == 2 %hold parameter 2 constant
                estimatedParam1 = fminsearch(@(param1)erfLikelihood([param1 param2], dataInRange), estimatedParams(1), options);
                estimatedParams = [estimatedParam1 param2];
            else
                estimatedParams = fminsearch(@(nlParams)erfLikelihood(nlParams, dataInRange), estimatedParams, options);
            end
        end
        
        if iter > 15
            disp(sprintf('robust fitting shouldn''t require %0.0f iterations -- looping?',iter)) ;        
        end
        
        % re-determine which data points lie in range based on new fit
        projectedProbs = 0.5 + 0.5*erf(estimatedParams(1)*data(1,:)+estimatedParams(2));
        
        aboveMin = projectedProbs > lowerProbLim;
        belowMax = projectedProbs < upperProbLim;
        
        minMovInRangeInd = find(aboveMin,1);
        maxMovInRangeInd = find(belowMax,1,'last');
        
        firstMovIndNew = max([minMovInRangeInd-1, 1]);
        lastMovIndNew = min([maxMovInRangeInd+1, length(projectedProbs)]);
        
        if (firstMovIndNew == firstMovInd && lastMovIndNew == lastMovInd)
            exitRobustFitting = true;
        else
            firstMovInd = firstMovIndNew;
            lastMovInd = lastMovIndNew;
        end
    end
end

a = estimatedParams(1);
b = estimatedParams(2);


xProj = min(data(1,:)):0.001:max(data(1,:));
projection = 0.5 + 0.5*erf(a*xProj+b);
projectionComplete = [xProj; projection]; 

error = norm(0.5 + 0.5*erf(a*data(1,:)+b) - data(2,:));


markerSize = 30*data(3,:)/max(data(3,:));

if makePlot
    figure('position', [400 400 600 400])
    hold on
    if any(data(1,:) < 0) %working in log scale
        title(['error = ', num2str(error), 10, 'threshold = ', num2str(exp(-b/a))])
    else %probably in linear scale
        title(['error = ', num2str(error), 10, 'threshold = ', num2str(-b/a) ', SD = ' num2str(1/(a*sqrt(2)))])
    end
    if ~isempty(lockedAmps)
        for i = 1:size(data, 2)
            if lockedAmps(i)
                plot(data(1,i), data(2,i),'k.', 'markerSize', markerSize(i));
            else
                plot(data(1,i), data(2,i),'.', 'markerSize', markerSize(i), 'markerFaceColor', [1 1 1], 'markerEdgeColor', [0.5 0.5 0.5]);
            end
        end
    else
        for i = 1:size(data, 2)
            plot(data(1,i), data(2,i),'k.', 'markerSize', markerSize(i), 'markerFaceColor', [1 1 1]);
        end
    end
    plot(xProj, projection,'m-');
    set(gca, 'ylim', [0 1])
    xlabel('current amplitude (µA)')
    ylabel('response rate')
    hold off
    
end

