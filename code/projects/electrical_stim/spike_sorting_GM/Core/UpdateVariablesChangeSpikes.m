function Gibbs = UpdateVariablesChangeSpikes(Gibbs,latenciesNew,neuronIndex)


Kn = Gibbs.params.Kn;
nNeurons = Gibbs.params.nNeurons;
J        = Gibbs.params.J;
I        = Gibbs.params.I;
Tfind0   = Gibbs.params.Tfind0;
E        = Gibbs.params.E;
dataVecJ = Gibbs.params.dataVecJ;
ActionPotentials = Gibbs.variables.ActionPotentials;


for n=neuronIndex
   
    for j = 1:J
        for i = 1:I(j)
           
            indLatency                 = find(Tfind0 == latenciesNew{n}(j,i));
            ActionPotentials{n,j}(i,:) = (Kn{n}(:,indLatency))';
        end
    end
end

[ResSpikesE ResSpikesEJ]=substractActionPotentials(dataVecJ,ActionPotentials,E);
  


Gibbs.variables.ActionPotentials = ActionPotentials;
Gibbs.variables.ResSpikesEJ      = ResSpikesEJ;
Gibbs.variables.ResSpikesE       = ResSpikesE;
  
     
Artifact = [];
for e = 1:E
        
        Gibbs    = sampleArtifact(Gibbs,e);
        Artifact = [Artifact Gibbs.variables.ArtifactE{e}];    
end

Gibbs.variables.Artifact = Artifact;


for e = 1:E
    for j=1:J
     
    Residuals{j,e}   = Gibbs.variables.ResSpikesEJ{e}{j}-repmat(Gibbs.variables.ArtifactE{e}(j,:),I(j),1);
        
    end
end
    

Gibbs.variables.Residuals = Residuals;
Gibbs                     = samplesigma(Gibbs);
Gibbs                     = LogisticRegression(Gibbs);
