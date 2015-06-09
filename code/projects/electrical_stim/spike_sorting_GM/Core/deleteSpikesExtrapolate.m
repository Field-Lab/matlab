function Gibbs = deleteSpikesExtrapolate(Gibbs,CondsDel)



spikes           = Gibbs.variables.spikes;
latencies        = Gibbs.variables.latencies;
ActionPotentials = Gibbs.variables.ActionPotentials;

nNeurons  = Gibbs.params.nNeurons;
I         = Gibbs.params.I;
J         = Gibbs.params.J;
T         = Gibbs.params.T;
E         = Gibbs.params.E;
dataVecJ  = Gibbs.params.dataVecJ;
Tfind0    = Gibbs.params.Tfind0;

for n = 1:nNeurons
    
    for j = CondsDel{n}
            spikes{n}(j,1:I(j))     = 0;
            latencies{n}(j,1:I(j))  = 0;
            ActionPotentials{n,j}   = zeros(I(j),E*T);
    end
end



Gibbs.variables.spikes    = spikes;
Gibbs.variables.latencies = latencies;


  
[ResSpikesE ResSpikesEJ]  = substractActionPotentials(dataVecJ,ActionPotentials,E);  
  

Gibbs.variables.ActionPotentials = ActionPotentials;
Gibbs.variables.ResSpikesEJ      = ResSpikesEJ;
Gibbs.variables.ResSpikesE       = ResSpikesE;
  
     
Artifact=[];
for e=1:E
        Gibbs    = SampleConditionalArtifact(Gibbs,e,CondsDel{e},setdiff([1:max(CondsDel{e})],CondsDel{e}));   
        Artifact = [Artifact Gibbs.variables.ArtifactE{e}];
   
end

Gibbs.variables.Artifact = Artifact;

% 
% for e = 1:E
%     for j=1:J
%      
%     Residuals{j,e}   = Gibbs.variables.ResSpikesEJ{e}{j}-repmat(Gibbs.variables.ArtifactE{e}(j,:),I(j),1);
%         
%     end
% end
    

Gibbs.variables.Residuals = Residuals;
%Gibbs                     = samplesigma(Gibbs);
Gibbs                     = LogisticRegression(Gibbs);
Gibbs                     = GibbsSamplerSpikesArtifact(Gibbs);
