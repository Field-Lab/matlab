function input = fillDefaultValues(input)
%  fillDefaultValues  sets the default values for input structure (and subsequent stages of algorithm).
%  it has to be executed after loadData for practical considerations, All changes in this default values
%  should be made after the execution of fillDefaultValues. Also, it sets all relevant variables after loadData (this means, 
%  this function could be somehow merged with laodData in order to leave it exclusively for setting of the default values)
%   Input:    - input: input structure
%            
%   Output:   - input: input structure with default values and remaining relevant
%               variables. See inline comments for details about those default values.
%               
%    
% Gonzalo Mena 6/2015 




% Set up defaults for optional parameters%

input.tracesInfo.J = size(input.tracesInfo.data,1); % Number of conditions (stimulus amps)
input.tracesInfo.E = length(input.tracesInfo.recElecs); % Number of recording electrodes
input.tracesInfo.T = size(input.tracesInfo.data{1,1},2); % Number of samples

for j=1:input.tracesInfo.J 
    I(j) = size(input.tracesInfo.data{j,1},1); % For each amplitude condition, the number of trials may or may not be the same. Therefore the number of trials is a function of J
end
input.tracesInfo.I = I;
input.neuronInfo.nNeurons = length(input.neuronInfo.neuronIds);



%Defaul sample rate

input.params.sampRate = 20000;

% Defines how many conditions must be within a particular breakpoint to use
% an artifact model of various degrees. length(degPolRule) = max polynomial
% fit minus one. For example, [1 50] uses a 0 order model for all cases
% between breakpoints with <=1 conditions, and a 1 order (linear) model for
% all cases with between 2 and 50 conditions.
input.params.initial.degPolRule   =  [1 4 15 50]; % Default to include a cubic polynomial model to account 
                                        % for axon bundle activation. Future defaults might not need a 
                                        % cubic model if there is a better way to approximate the axon 
                                        % bundle activation shape.%[1 50]; 
input.params.initial.Tdivide      =  [0 floor(input.tracesInfo.T/2) input.tracesInfo.T]; % After the initialization, a regularized
                                                                                         % model for the artifact is built. Tdivide divides the time window
                                                                                         % we are looking at in different chunks, 
                                                                                         % and each chunk is treated separately, as the artifact features change over time 
                                                                                         %(at the begining oscilations are wilder, for example)
                                                                                         % In this case, we are dividing the observation window in two chunks to find separate
                                                                                         % time regularization hyperparameters for each.

                                                                                         
                                                                                         
 
                                                                                         
  
%Default minimum and maximum recording times. If the first part of the
%traces are not relevant, it could be worth trying make Tmin>1

Tmin     = input.tracesInfo.Trange(1);
Tmax     = input.tracesInfo.Trange(2);
 
                                                                                                                                                                         
input.params.initial.TfindRange        =    [Tmin+5 Tmax-5]; % Range of times for which we actually search for spikes
input.params.initial.TfindRangeRel     =    input.params.initial.TfindRange - input.tracesInfo.Trange(1) +1; %Range of times to look for spikes, relative to the first                                                                                                          %time first sample being considered
input.params.initial.rho               =    0; % Initialization parameter. Degree of sparseness that is imposed on the initial solution.

% parameters imposed on the Newton method for finding the regularization hyperparameters lambda
input.params.initial.Newton.beta  = 0.5; %beta and a correspond to the gradient descent method with backtracking linesearches
input.params.initial.Newton.a     = 0.5;
input.params.initial.Newton.thres = 0.01; % Iterations will stop once the gradient norm is smaller than this value.



% don't change, decreases the sensitivity of the sigmoidal fit to low or
% high probability spiking (low probaiblity of spiking is defined by this value)
% This means, for computation of activation curves, very rare spike is ommitted
% and if spiking is very close to probability one then it is treated as actually having probability one
input.params.Gibbs.thresLogistic = 0.15*ones(1,input.tracesInfo.J);


%Maximum number of Gibbs sampling iterations. This number is important, in
%some cases we have to set it small so we won't waste time iterating and
%getting the same results. The problem is that if it is too small then the
%Gibbs sampler won't 'move' enough
input.params.maxIterGibbs = 10;

% heuristic parameters - addresses cases of high/low activation in
% subsequent stimulation conditions.
input.params.Heuristic.ResampleAlltimes    = 1;    %If 1, then intrapolation or resampling of the artifact is done taking in consideration the entire recording window. Otherwise resampling will be done at the pieces defined by Tdivide
input.params.Heuristic.ExtrapolateAlltimes = 1; %If 1, then extrapolation of the artifact is done taking in consideration the entire recording window. Otherwise resampling will be done at the pieces defined by Tdivide
input.params.Heuristic.ActivationThres  = 0.3;  %Definition of the activation threshold 
input.params.Heuristic.VarianceFactor   = 1.2;  %The value such that if deletion increases the standard deviation by this factor at some condition, deletion will stop
input.params.Heuristic.LogisticRegPValueThres = 10^-4;  %For assesing lack of fit in logistic regression a significance level  has to be specified. If the chi square test results in a value smaller than this significance level, a lack of fit is detected 
input.params.Heuristic.maxCondsChange         = 3;   %Used by HExtrapolateBadLogRegression and HResampleBadLogRegression,  maximum number of simultaneous 'worst fit' conditions that will be taken into account for resampling or extrapolating.
input.params.Heuristic.numErrorConds          = 4;   %Used by HExtrapolateBadLogRegression and HResampleBadLogRegression,  number of 'worst fit' conditions that will be taken into account for resampling or extrapolating.
input.params.Heuristic.maxCondsShift          = 2;   %Used by HResampleIfLackOfSpiking. It is the maximum number of conditions that will be exprapolated simultaneously
input.params.Heuristic.maxIter                = 50;   %Maximum number of iterations of the heuristic HExtrapolateBadLogRegression
input.params.Heuristic.HighThres              = 0.8;  %Definition of a high spike probability


input.params.Heuristic.LowThres               = 0.25; %Definition of a low spike probability
input.params.Heuristic.VeryLowThres           = 0.1; %Definition of a very low spike probability



