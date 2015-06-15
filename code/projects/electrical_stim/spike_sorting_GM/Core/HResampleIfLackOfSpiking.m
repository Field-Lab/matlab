function [Gibbs Log] =  HResampleIfLackOfSpiking(Gibbs,input,Log,neuronIndex,varargin)


ActivationThres        = input.params.Heuristic.ActivationThres;
maxShift               = input.params.Heuristic.maxCondsChange;
LowThres               = input.params.Heuristic.LowThres;

E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
I           = Gibbs.params.I;
nNeurons    = Gibbs.params.nNeurons;
breakPoints = input.tracesInfo.breakPoints;
nActiveNeurons = length(neuronIndex);
prefElectrodes = input.neuronInfo.prefElectrodes;


for n=neuronIndex
    contLog(n) = Log(n).params.contLogHeuristic;
end

if(nargin>=5)
    CondsInit=varargin{1};
end

if(nargin == 6)
    maxCond = varargin{2};
else
    maxCond = J;
end

for e=1:E
    CondsResampleOld{e} =[];
    CondsResample{e}    =[];
    shiftRight(e)=0;
end



for n=1:nActiveNeurons
    aNeuron = neuronIndex(n);
    CondsResample{n} = [];
    
    if(nargin>=5)
        e = prefElectrodes{aNeuron}(1);
        if(~isempty(find((CondsInit{aNeuron}-1)==breakPoints{e})))
            CondsInit{aNeuron}=CondsInit{aNeuron}+1;
        end
    end

end


contIter=1;
while(true)
    clear CondsAux
    for n=1:nNeurons
        CondsAux{n}=[];
    end
    
    spikes = Gibbs.variables.spikes;
    flags = ones(1,nActiveNeurons);
   
    for n=1:nActiveNeurons
        aNeuron = neuronIndex(n);
        
        auxLow  = find((nansum(spikes{aNeuron}')./I)<=ActivationThres);
        auxHigh = find((nansum(spikes{aNeuron}')./I)>ActivationThres);
        
        if(nargin>=5)
            if(isempty(auxHigh)||contIter==1)
                CondsAux{aNeuron} = CondsInit{aNeuron};
                flags(n)          = isempty(CondsAux{aNeuron});
                continue
            end
        end
        if(nargin<5||contIter>1)
            if(isempty(auxHigh)) 
                continue
            else
                auxLow=auxLow(auxLow>auxHigh(end));
                if(isempty(auxLow))
                    
                    continue;
                else
                    CondsAux{aNeuron}=auxLow(1);
                    flags(n) = 0;
                end
            
            end
        end       
    end
    
        
    
    if(prod(flags)==1)
        return;
    end
    

    CondsAux    = CondsNeuron2Electrode(input,CondsAux);
    
    flags = zeros(1,E);
    for e=1:E
        if(isempty(CondsAux{e}))
            CondsResample{e}=[];
            flags(e)=1; 
            continue
        else
            CondsResample{e}=[CondsAux{e}:CondsAux{e}+shiftRight(e)];
            if(isequal(CondsResample{e},CondsResampleOld{e}))
                shiftRight(e)=shiftRight(e)+1;
                CondsResample{e}=[CondsAux{e}:CondsAux{e}+shiftRight(e)];
            else
                shiftRight(e)=0;
            end
            if(shiftRight(e)>maxShift)
                flags(e)=1;
                shiftRight(e)    = maxShift;
     
            end
            
            if(~isempty(find((breakPoints{e}+1)==CondsResample{e}(1))));
                CondsResample{e}=[];
                flags(e)=1;
                continue
            end
            
            if(CondsResample{e}(end)>maxCond)
                CondsResample{e}=[];
                flags(e)=1;
                
                continue
            end

        end
        
    end
    
    
    if(prod(flags)==1)
        return;
    end
    
    CondsResampleOld = CondsResample;
    
    for n=1:neuronIndex
        e = prefElectrodes{n}(1);
        if(flags(e) == 0)
            contLog(n)=contLog(n)+1;
            Log(n).Heuristic{contLog(n)}=['Going to extrapolate at conditions' num2str(CondsResample{e}) ' since there is low spiking following activation'];
            Log(n).params.contLogHeuristic = contLog(n);     
        end
    end


        
  
    Artifact=[];
    
    for e=1:E
        
        if(flags(e)==0)
            Gibbs       = SampleConditionalArtifact(Gibbs,e,CondsResample{e},setdiff([1:max(CondsResample{e})],CondsResample{e})); 
        end
        Artifact    = [Artifact Gibbs.variables.ArtifactE{e}];
    end
    
    Gibbs.variables.Artifact = Artifact;
    Gibbs                    = GibbsSamplerSpikesArtifact(Gibbs);
    
    contIter=contIter+1;
end
    