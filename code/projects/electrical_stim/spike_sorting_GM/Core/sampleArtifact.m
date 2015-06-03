function Gibbs = sampleArtifact(Gibbs,e)



Lambdas     = Gibbs.params.Lambdas;
Xj          = Gibbs.params.Xj;
X           = Gibbs.params.X;    
I           = Gibbs.params.I;
J           = Gibbs.params.J;
T           = Gibbs.params.T;

Beta          = Gibbs.variables.Beta;
sigma         = Gibbs.variables.sigma;
ResSpikesEJ   = Gibbs.variables.ResSpikesEJ;      
ResSpikesE    = Gibbs.variables.ResSpikesE;       




SigmaInv=[];

for j=1:J
    SigmaInv = sparse(blkdiag(SigmaInv,1/sigma(e,j)^2*speye(I(j)*T)));
end

SigmaPost{e}       =    (Lambdas{e}+X'*SigmaInv*X)^(-1);
SigmaPost{e}       =    (SigmaPost{e}+SigmaPost{e}')/2;
muBetaPost(:,e)    =    SigmaPost{e}*X'*SigmaInv*ResSpikesE{e};


Beta(:,e)          =    mvnrnd(muBetaPost(:,e)',full(SigmaPost{e}))';

for j=1:J
    Gibbs.variables.ArtifactE{e}(j,:) =  (Xj{j}*Beta(:,e))';
end

Gibbs.variables.Beta         =  Beta;
Gibbs.variables.SigmaPost{e} =  (SigmaPost{e}+SigmaPost{e}')/2;


