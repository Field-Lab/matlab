function [Gibbs Log]=HResampleBadLogRegression(input,Gibbs,Log)


type        = input.params.Heuristic.ResampleAlltimes;
E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
nNeurons    = Gibbs.params.nNeurons;

thres          = input.params.Heuristic.LogisticRegPValueThres;
maxCondsChange = input.params.Heuristic.maxCondsChange;



flags = zeros(nNeurons,1);

for n=1:nNeurons
    contLog(n) = Log(n).params.contLogHeuristic;
    CondsChangeold{n}=[];
    nCondsChange(n)=1;
end

while(prod(flags)==0)

    p2          = Gibbs.diagnostics.Logisticsp2GoF;
    condsSorted = Gibbs.diagnostics.condSortError;
        
    for n=1:nNeurons

        if(p2(n)>thres)
            
            CondsChange{n}=[];
            CondsChangeold{n}=CondsChange{n};
            flags(n)=1;
            continue
    
        end


        CondsChange{n}=condsSorted(n,1:min(nCondsChange(n),maxCondsChange));

        if(length(CondsChange{n})==maxCondsChange)
            random          = unidrnd(maxCondsChange-1);
            nCondsChange(n) = random;
            randSample      = sort(randsample(maxCondsChange,random))';
            CondsChange{n}  = condsSorted(n,randSample);
            flags(n)=1;
    
            continue
        elseif(isequal(unique(CondsChange{n}),unique(CondsChangeold{n})))
            nCondsChange(n)=nCondsChange(n)+1;
            CondsChange{n}=condsSorted(n,1:min(nCondsChange(n),maxCondsChange));

        else
            nCondsChange(n)=1;
            CondsChange{n}=condsSorted(n,1:min(nCondsChange(n),maxCondsChange));

        end


    CondsChangeold{n}=CondsChange{n};

    end

    CondsChange    = CondsNeuron2Electrode(input,CondsChange);
    
    for n=1:nNeurons
        if(~isempty(CondsChange{n}))
            contLog(n)=contLog(n)+1;
            Log(n).Heuristic{contLog(n)}=['Logistic regression fit is poor at Conditions ' num2str(CondsChange{n}) ' p = ' num2str(Gibbs.diagnostics.Logisticsp2GoF(n))];
            Log(n).params.contLogHeuristic = contLog(n);
        end
    end
    ArtifactNew=[];

    for e = 1:E
    
        if(type ==1)
            
            Gibbs       = SampleConditionalArtifact(Gibbs,e,CondsChange{e},setdiff([1:J],CondsChange{e}));
        else
            Gibbs       = SampleConditionalArtifactT(Gibbs,e,CondsChange{e},setdiff([1:J],CondsChange{e})); 
        end
        ArtifactNew = [ArtifactNew Gibbs.variables.ArtifactE{e}];
   
   
    end
    Artifact                 = ArtifactNew;
    Gibbs.variables.Artifact = ArtifactNew;
    Gibbs                    = GibbsSamplerSpikesArtifact(Gibbs);

end
