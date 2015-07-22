function Gibbs =SampleConditionalArtifactT(Gibbs,e,CondsChange,CondsOn)


T           =  Gibbs.params.T;
Tdivide     =  Gibbs.params.Tdivide;
J           =  Gibbs.params.J;
LambdaInv   =  Gibbs.params.LambdasInv{e};
Xj          =  Gibbs.params.Xj;
Beta        =  Gibbs.variables.Beta;
sigmaT      =  Gibbs.variables.sigmaT;

if(isempty(CondsChange))
  return
end


for j = 1:length(CondsChange)
    
    clear aux
    for t = 1:length(Tdivide)-1
           
        sigmat(t) = sigmaT{t}(e,CondsChange(j)).^2;
             
    end
    sigmatprop = sigmat/sum(sigmat);
    indT(j)    = find(mnrnd(1,sigmatprop,1)==1);
          
end


mu = zeros(size(Beta,1),1);

indsBetaChange=[];
CondsOn=setdiff(CondsOn,CondsChange);

for j = 1:length(CondsChange)
 
    indsTdivide     = [(Tdivide(indT(j))+1):Tdivide(indT(j)+1)];
    indsBetaChange = [indsBetaChange (CondsChange(j)-1)*T+indsTdivide];

end

indsBetaOn = [];
for j = 1:length(CondsOn)
    indsBetaOn = [indsBetaOn [1:T]+(CondsOn(j)-1)*T];
end

for j = 1:length(CondsChange)

      indsTdivide     = [(Tdivide(indT(j))+1):Tdivide(indT(j)+1)];
      indsTdivideNot  = setdiff([1:T],indsTdivide);
      indsBetaOn = [indsBetaOn (CondsChange(j)-1)*T+indsTdivideNot];

end

indsBetaChange = sort(indsBetaChange);
indsBetaOn     = sort(indsBetaOn);
indsUnion      = union(indsBetaChange,indsBetaOn);

LambdaInv = LambdaInv(indsUnion,indsUnion);
mu        = mu(indsUnion);
Betanew   = Beta(:,e);
x         = Betanew(indsBetaOn);


CondsUnion = union(CondsChange,CondsOn);

for j = 1:length(CondsChange)
    CondsChangeRel(j) = find(CondsChange(j)==CondsUnion);
end

indsBetaChangeRel=[];

for j = 1:length(CondsChangeRel)

    indsTdivide=[(Tdivide(indT(j))+1):Tdivide(indT(j)+1)];
    indsBetaChangeRel = [indsBetaChangeRel (CondsChangeRel(j)-1)*T+indsTdivide];

end


[BetaCond] = SampleConditionalNormal(mu,LambdaInv,x,indsBetaChangeRel);

Betanew(indsBetaChange)=BetaCond;


for j=1:J;
    
    ArtifactE{e}(j,:) = (Xj{j}*Betanew)';
    
end
Gibbs.variables.Beta(:,e)    = Betanew;
Gibbs.variables.ArtifactE{e} = ArtifactE{e};

