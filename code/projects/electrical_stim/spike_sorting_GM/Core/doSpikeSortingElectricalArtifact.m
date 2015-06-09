function [Gibbs GibbsNoDelete initial input Log] = doSpikeSortingElectricalArtifact(input,initial)


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

Gibbs = GibbsSamplerSpikesArtifact(Gibbs);

[Gibbs Log] = deleteSpikesBeginning(input,Gibbs,Log,1);
[Gibbs Log] = HResampleBadLogRegression(input,Gibbs,Log);
[Gibbs Log] = HResampleIfLackOfSpiking(Gibbs,input,Log,[1:nNeurons]);
[Gibbs Log] = HResampleAboveActivation(Gibbs,input,Log);
[Gibbs Log] = deleteSpikesBeginning(input,Gibbs,Log,1);

GibbsNoDelete = Gibbs;

[Gibbs Log] = HDeleteSpikesAndResample(GibbsNoDelete,input,Log);
[Gibbs Log] = deleteSpikesBeginning(input,Gibbs,Log,1);
[Gibbs Log] = HaddSpikesResample(Gibbs,input,Log);
[Gibbs Log] = HResampleAboveActivation(Gibbs,input,Log);


end




