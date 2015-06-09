function Gibbs = GibbsSamplerSpikesArtifact(Gibbs)
    

dataVecJ    = Gibbs.params.dataVecJ;
TfindRel0   = Gibbs.params.TfindRel0;
Tfind0      = Gibbs.params.Tfind0;
lambdas     = Gibbs.params.lambdas;
Xj          = Gibbs.params.Xj;
X           = Gibbs.params.X;    
Kn          = Gibbs.params.Kn;   
I           = Gibbs.params.I;
E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
nNeurons    = Gibbs.params.nNeurons;


sigma            = Gibbs.variables.sigma;
Artifact         = Gibbs.variables.Artifact; 
Probs            = Gibbs.variables.Probs;     
ActionPotentials = Gibbs.variables.ActionPotentials;

contIter = 1;

while(true)
  
    if(contIter>1)
        spikesold = Gibbs.variables.spikes;
    end
    
    Gibbs = sampleSpikes(Gibbs);

    Artifact=[];
    
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
   
    Gibbs = samplesigma(Gibbs);
    Gibbs = LogisticRegression(Gibbs);


    if(contIter>1)
        changeSpikes=0;
        for n=1:nNeurons
            for j = 1:J
            changeSpikes=changeSpikes+sum(abs(Gibbs.variables.spikes{n}(j,1:I(j))-spikesold{n}(j,1:I(j))));
        
            end
        end
        
        if(changeSpikes==0)
            return;
        end
    end
contIter = contIter+1;
end
