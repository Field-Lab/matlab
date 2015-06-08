function [Gibbs]=SampleConditionalArtifact(Gibbs,e,CondsChange,CondsOn)

T           =  Gibbs.params.T;
J           =  Gibbs.params.J;
LambdaInv   =  Gibbs.params.LambdasInv{e};
Xj          =  Gibbs.params.Xj;
Beta        =  Gibbs.variables.Beta;

if(isempty(CondsChange))
  return
end
    
mu = zeros(size(Beta,1),1);

indsBeta1 = [];
for j=1:length(CondsOn)
    indsBeta1 = [indsBeta1 [1:T]+(CondsOn(j)-1)*T];
end

indsBeta2 = [];

for j = 1:length(CondsChange)

    indsBeta2 = [indsBeta2 [1:T]+(CondsChange(j)-1)*T];
end


indsBeta    = union(indsBeta1,indsBeta2);
LambdaInv   = LambdaInv(indsBeta,indsBeta);

mu = mu(indsBeta);

Betanew = Beta(:,e);

x = Betanew(indsBeta1);
condsUnion = union(CondsChange,CondsOn);

for j = 1:length(CondsChange)
    CondsChangeRel(j) = find(CondsChange(j)==condsUnion);
end

indsBeta2Rel = [];

for j = 1:length(CondsChangeRel)
    
    indsBeta2Rel = [indsBeta2Rel [1:T]+(CondsChangeRel(j)-1)*T];
end


[BetaCond] = SampleConditionalNormal(mu,LambdaInv,x,indsBeta2Rel);


Betanew(indsBeta2)=BetaCond;


for j=1:J;
    
    ArtifactE{e}(j,:) = (Xj{j}*Betanew)';
    
end
Gibbs.variables.Beta(:,e)    = Betanew;
Gibbs.variables.ArtifactE{e} = ArtifactE{e};

