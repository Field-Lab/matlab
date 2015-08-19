function Gibbs = extrapolateFromCondition(Gibbs,input,neuronIndex,CondInit,breakRange)



E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
I           = Gibbs.params.I;
nNeurons    = Gibbs.params.nNeurons;
prefElectrodes = input.neuronInfo.prefElectrodes;
e              = prefElectrodes{neuronIndex}(1);
breakPoints    = [input.tracesInfo.breakPoints{e} J];
LastCond       = breakPoints(breakRange);


Cond           = CondInit;

while(Cond <= LastCond)
    
    for n=1:nNeurons
        Conds{n}=[];
    end
    Conds{neuronIndex}=Cond;
    Conds    = CondsNeuron2Electrode(input,Conds);
    
    Artifact=[];
    
    for e=1:E
        
        if(~isempty(Conds{e}))
            Gibbs       = SampleConditionalArtifact(Gibbs,e,Conds{e},setdiff([1:max(Conds{e})],Conds{e})); 
        end
        Artifact    = [Artifact Gibbs.variables.ArtifactE{e}];
    end
    
    Gibbs.variables.Artifact = Artifact;
    Gibbs                    = GibbsSamplerSpikesArtifact(Gibbs);
    Cond=Cond+1;
end