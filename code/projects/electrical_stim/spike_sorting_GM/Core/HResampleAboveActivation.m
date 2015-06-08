function [Gibbs Log] = HResampleAboveActivation(Gibbs,input,Log)

ActivationThres = input.params.Heuristic.ActivationThres;
HighThres       = input.params.Heuristic.HighThres;
E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
I           = Gibbs.params.I;
nNeurons    = Gibbs.params.nNeurons;
prefElectrodes = input.neuronInfo.prefElectrodes;
breakPoints    = input.tracesInfo.breakPoints;

for e=1:E
    CondsResampleOld{e}=[];
end

for n=1:nNeurons
    contLog(n) = Log(n).params.contLogHeuristic;
    CondsResample{n} = [];
end

while(true)
    
    spikes = Gibbs.variables.spikes;
    CondsActivation = detectActivation(input,spikes,ActivationThres);
    for n=1:nNeurons
        
        aux=find((nansum(spikes{n}')./I)>HighThres);
        if(~isempty(aux))
            if(aux(1)<=CondsActivation(n))
                CondsResample{n}=aux(1);
            end
        end
    end
    
    CondsResample    = CondsNeuron2Electrode(input,CondsResample);
    
    flags = zeros(E,1);
    
    for e = 1:E
        if(isempty(CondsResample{e}))
            flags(e) = 1;
        else
            if(isequal(CondsResample{e},CondsResampleOld{e})||~isempty(find(CondsResample{e}==(breakPoints{e}+1))))
                flags(e)=1;
            end
        end
    end
    
    if(prod(flags)==1)
        return;
    end
    for n = 1:nNeurons
        e = prefElectrodes{n}(1);
        if(flags(e) == 0)
            contLog(n)=contLog(n)+1;
            Log(n).Heuristic{contLog(n)}=['Condition ' num2str(CondsResample{e}) ' has a lot activation. Going to extrapolate at that condition'];
            Log(n).params.contLogHeuristic = contLog(n);
        end
    end
    
    
    CondsResampleOld = CondsResample;
    
    
    Artifact=[];
    
    for e=1:E
        
        Gibbs       = SampleConditionalArtifact(Gibbs,e,CondsResample{e},setdiff([1:max(CondsResample{e})],CondsResample{e}));
        Artifact    = [Artifact Gibbs.variables.ArtifactE{e}];
    end
    
    Gibbs.variables.Artifact = Artifact;
    Gibbs                    = GibbsSamplerSpikesArtifact(Gibbs);
    
    
end