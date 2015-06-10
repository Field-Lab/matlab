function input = fillDefaultValues(input)


input.tracesInfo.J = size(input.tracesInfo.data,1); % Number of conditions (stimulus amps)
input.tracesInfo.E = length(input.tracesInfo.recElecs); % Number of recording electrodes
input.tracesInfo.T = size(input.tracesInfo.data{1,1},2); % Number of samples

for j=1:input.tracesInfo.J 
    I(j) = size(input.tracesInfo.data{j,1},1); % For each amplitude condition, the number of trials may or may not be the same. Therefore the number of trials is a function of J
end
input.tracesInfo.I = I;
input.neuronInfo.nNeurons = length(input.neuronInfo.neuronIds);
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
input.params.initial.Tdivide      =  [0 floor(input.tracesInfo.T/2) input.tracesInfo.T];


Trange = input.tracesInfo.Trange;
input.params.initial.TfindRange        =    [Trange(1)+5 Trange(2)-5]; % Range of times for which we actually search for spikes
input.params.initial.TfindRangeRel     =    input.params.initial.TfindRange - input.tracesInfo.Trange(1) +1;
input.params.initial.rho               =    0; % Initialization parameter. Degree of sparseness that is imposed on the initial solution.

% parameters imposed on the initial solution
input.params.initial.Newton.beta  = 0.5;
input.params.initial.Newton.a     = 0.5;
input.params.initial.Newton.thres = 0.01;

% don't change, decreases the sensitivity of the sigmoidal fit to low
% probability spiking. 
input.params.Gibbs.thresLogistic = 0.15*ones(1,input.tracesInfo.J);

% heuristic parameters - addresses cases of high/low activation in
% subsequent stimulation conditions.
input.params.Heuristic.ResampleAlltimes = 1;
input.params.Heuristic.ExtrapolateAlltimes = 1;
input.params.Heuristic.ActivationThres  = 0.3;
input.params.Heuristic.VarianceFactor   = 1.2;
input.params.Heuristic.LogisticRegPValueThres = 10^-4;
input.params.Heuristic.maxCondsChange         = 3;
input.params.Heuristic.numErrorConds          = 4;
input.params.Heuristic.maxCondsShift          = 2;
input.params.Heuristic.maxIter                = 50;
input.params.Heuristic.HighThres              = 0.8;
input.params.Heuristic.LowThres               = 0.15;
input.params.Heuristic.VeryLowThres           = 0.1;



