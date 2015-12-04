function [Gibbs GibbsNoDelete initial input Log] = doSpikeSortingElectricalArtifact(input,initial)
%doSpikeSortingElectricalArtifact takes the input and initial structures as input and with that, executes the core
% part of the spike sorting algorithm: In the first place, it does Gibbs sampling iterations (iterating between
% spikes-latencies, Artifact, residual variances and logistic regression for spikes)
% until no changes are observed, or until a maximum number of iterations,
% specified in input.params.Gibbs.maxIter is exceeded (the choice of this number can be
% critical for the several neurons-electrodes context. Variables are stored in the Gibbs structure, and they suffer changes from
% iteration to iteration. After these first Gibbs sampler iterations, doSpikeSortingElectricalArtifact
% will apply sequentially the different heuristics (files that start with a H) in case some problem has been diagnosed
% (for example, lack of fit with the underlying logistic regression model)
% Almost all of these heuristics are based in re-estimating (extrapolating or interpolating)
% the artifact in conditions for which problems have been diagnosed, and doing more Gibbs iterations again, until
% no observed changes in spikes are observed (or until iterations exceed the maximum allowed).
% Nonetheless, one of the heuristics is more aggressive  and will try to delete ALL spikes in certain conditions, until this 
% delation leads to an increase in the residual variance. This strategy seems to work well in practice,
%and in the forst case, the solution withouth deletion will be part of the output as well.
% Notice that sequentiallness of Heuristics was done here only because it seemed that was the right order
%However, there are possibilities in exploring how the different Heuristics should succedd the others
%or even, if they could be included in a while loop that will be active until no more changes are needed.
%But for now, keeping things like this should provide reasonable results
%
%Input:        input and initial structure
%Output:       -same input and initial structure (this could be changed with no harm)
%              -Gibbs: a structure such that Gibbs.variables contains all relevant variables that are sampled 
%                      Also, Gibbs.params contains the paramters for the Gibbs sampler (as maximum number of iterations)
%                      And Gibbs.diagnostics contain information about logistic regression fit to spikes (not relevant now)
%              -GibbsNoDelete: same as Gibbs, but results reflect what happened before iterations
%              -Log: a structure such that Log(n) is the Log for neuron with index n, and
%               Log(n).Heuristic is a array of strings containing information about what changes were by the heuristics
%               (including interpolations/extrapolations) of artifact and deletion
%               Log(n).Deletion is a vector containing information about if it was or no deletion for neuron n at the different breakpoints ranges.
%Important:     Essentially, spike sorting solutions correspond to Gibbs.variables.spikes and Gibbs.variables.latencies (the same with GibbsNoDelete, in case
%               the experimenter thinks they can give useful information



%%
%load relevant variables and create Gibbs structure
data           = input.tracesInfo.data;
dataVecJ       = input.tracesInfo.dataVecJ;
templates      = input.neuronInfo.templates;
TfindRangeRel  = input.params.initial.TfindRangeRel;
TfindRel       = [TfindRangeRel(1):TfindRangeRel(2)];
TfindRel0      = [0 TfindRel]; 
neuronIds      = input.neuronInfo.neuronIds;

E = input.tracesInfo.E;
T = input.tracesInfo.T;
I = input.tracesInfo.I;
J = input.tracesInfo.J;

nNeurons = input.neuronInfo.nNeurons;

for n=1:nNeurons
    Log(n).params.contLogHeuristic = 0;
    for j=1:J
        Gibbs.variables.spikes{n}(j,1:I(j))=100000;
    end
end

lengthDataVec   = sum(I)*E*T;



K = makeToeplitz(templates,TfindRel0,T);

for n=1:nNeurons
    
  
    indn=zeros(1,nNeurons);
    indn(n)=1;
    
    for t=1:length(TfindRel0)
        
        indt       = zeros(1,length(TfindRel0));
        indt(t)    = 1;
        indnt      = kron(indn,indt);
        Kn{n}(:,t) = K*indnt';
    end
end


Gibbs.params.I = I;
Gibbs.params.E = E;
Gibbs.params.J = J;
Gibbs.params.T = T;
Gibbs.params.nNeurons = nNeurons;
Gibbs.params.maxIterGibbs = input.params.maxIterGibbs;

Gibbs.params.Kn                     = Kn;
Gibbs.params.X                      = initial.params.X;
Gibbs.params.Xj                     = initial.params.Xj;

for e=1:E
    Gibbs.params.matricesReg(e).Prods = initial.params.matricesReg(e).Prods;
end

Gibbs.params.lambdas                = initial.params.lambdas;
Gibbs.params.Lambdas                = initial.params.Lambdas;
Gibbs.params.LambdasInv             = initial.params.LambdasInv;
Gibbs.params.TfindRel0              = TfindRel0;
Gibbs.params.Tfind0                 = [0 input.params.initial.TfindRange(1):input.params.initial.TfindRange(2)];
Gibbs.params.dataVecJ               = dataVecJ;
Gibbs.params.a0                     = initial.params.a0;
Gibbs.params.b0                     = initial.params.b0;
Gibbs.params.lambdaLogReg           = initial.params.lambdaLogReg;
Gibbs.params.alphaLogReg            = initial.params.alphaLogReg;
Gibbs.params.Tdivide                = input.params.initial.Tdivide;
Gibbs.params.thresLogistic          = input.params.Gibbs.thresLogistic;

Gibbs.variables.sigma            = initial.sigma;
Gibbs.variables.Artifact         = initial.Artifact;
Gibbs.variables.Probs            = initial.Probs;
Gibbs.variables.ActionPotentials = initial.ActionPotentials;  
Gibbs.variables.Beta             = zeros(size(Gibbs.params.Lambdas{1},1),E);     
%%

%Do first Gibbs sampler iterations
Gibbs = GibbsSamplerSpikesArtifact(Gibbs);

%%
%Do heuristics
[Gibbs, Log] = HdeleteSpikesBeginning(input,Gibbs,Log,1);
[Gibbs, Log] = HResampleBadLogRegression(input,Gibbs,Log);
[Gibbs, Log] = HResampleIfLackOfSpiking(Gibbs,input,Log,[1:nNeurons]);
[Gibbs, Log] = HResampleAboveActivation(Gibbs,input,Log);
[Gibbs, Log] = HdeleteSpikesBeginning(input,Gibbs,Log,1);
%Save solution before executing Heuristics
GibbsNoDelete = Gibbs;

[Gibbs, Log] = HDeleteSpikesAndResample(GibbsNoDelete,input,Log);
[Gibbs, Log] = HdeleteSpikesBeginning(input,Gibbs,Log,1);
[Gibbs, Log] = HaddSpikesResample(Gibbs,input,Log);
[Gibbs, Log] = HResampleAboveActivation(Gibbs,input,Log);


end




