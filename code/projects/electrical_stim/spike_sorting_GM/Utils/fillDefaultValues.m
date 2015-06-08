function input = fillDefaultValues(input)


input.tracesInfo.J = size(input.tracesInfo.data,1);
input.tracesInfo.E = length(input.tracesInfo.recElecs);
input.tracesInfo.T = size(input.tracesInfo.data{1,1},2);

for j=1:input.tracesInfo.J 
    I(j) = size(input.tracesInfo.data{j,1},1);
end
input.tracesInfo.I = I;
input.neuronInfo.nNeurons = length(input.neuronInfo.neuronIds);

input.params.sampRate                =    20000;

input.params.initial.degPolRule   =  [1 50];
input.params.initial.Tdivide      =  [0 floor(input.tracesInfo.T/2) input.tracesInfo.T];


Trange = input.tracesInfo.Trange;
input.params.initial.TfindRange        =    [Trange(1)+5 Trange(2)-5];
input.params.initial.TfindRangeRel     =    input.params.initial.TfindRange - input.tracesInfo.Trange(1) +1;
input.params.initial.rho               =    0;

input.params.initial.Newton.beta  = 0.5;
input.params.initial.Newton.a     = 0.5;
input.params.initial.Newton.thres = 0.01;

input.params.Gibbs.thresLogistic = 0.15*ones(1,input.tracesInfo.J);



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


