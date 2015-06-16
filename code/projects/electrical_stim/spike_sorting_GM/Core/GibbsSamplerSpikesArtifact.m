function Gibbs = GibbsSamplerSpikesArtifact(Gibbs)
%GibbsSamplerSpikesArtifact  updates the Gibbs.variables field, by doing swipes
% between sampling spikes and latencies, Artifact, residual variances and logistic regression fit
% It stops until no more changes in Gibbs.variables.spikes are seen or until the maximum number of iterations
% is exceeded

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
maxIterGibbs = Gibbs.params.maxIterGibbs;

sigma            = Gibbs.variables.sigma;
Artifact         = Gibbs.variables.Artifact; 
Probs            = Gibbs.variables.Probs;     
ActionPotentials = Gibbs.variables.ActionPotentials;

contIter = 1;

while(contIter<=maxIterGibbs)
  
    spikesold = Gibbs.variables.spikes;
    
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
    
    flags=zeros(nNeurons,1);
    
        
        for n=1:nNeurons
            changeSpikes(n)=0;
            for j = 1:J
                changeSpikes(n)=changeSpikes(n)+sum(abs(Gibbs.variables.spikes{n}(j,1:I(j))-spikesold{n}(j,1:I(j))));
                
            end
            if(changeSpikes(n)<=1)
                flags(n)=1;
            end
            
            
        end
        if(prod(flags)==1)
            return;
        end
        contIter = contIter+1;
    end
