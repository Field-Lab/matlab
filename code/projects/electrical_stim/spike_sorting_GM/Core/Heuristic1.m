function Gibbs=Heuristic1(Gibbs)




E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
nNeurons    = Gibbs.params.nNeurons;

thres       = 10^-4;
maxCondsChange = 3;


flags = zeros(nNeurons,1);

for n=1:nNeurons
    CondsChangeold{n}=[];
    nCondsChange(n)=1;
end

while(prod(flags)==0)

    p1          = Gibbs.diagnostics.Logisticsp1GoF;
    condsSorted = Gibbs.diagnostics.condSortError;
        
    for n=1:nNeurons

        if(p1(n)>thres)
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

    ArtifactNew=[];

    for e= 1:E
    [e CondsChange{e}]
        Gibbs       = SampleConditionalArtifact(Gibbs,e,CondsChange{e},setdiff([1:J],CondsChange{e}));
        ArtifactNew = [ArtifactNew Gibbs.variables.ArtifactE{e}];
   
   
    end
    Artifact                 = ArtifactNew;
    Gibbs.variables.Artifact = ArtifactNew;
    Gibbs                    = GibbsSamplerSpikesArtifact(Gibbs);

end
