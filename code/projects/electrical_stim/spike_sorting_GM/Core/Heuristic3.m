function Gibbs = Heuristic3(Gibbs)
E           = Gibbs.params.E;
J           = Gibbs.params.J;
nNeurons    = Gibbs.params.nNeurons;

thres                  = 10^-4;
maxIter                = 1;
nIter                  = 0;



flags=zeros(nNeurons,1);
while(prod(flags)==0&&nIter<=maxIter)
    flags       = zeros(nNeurons,1);
    p1          = Gibbs.diagnostics.Logisticsp1GoF;
    
    if(mod(nIter,2)==0)
        Jindex = [1:J];
    else
        Jindex = -[-J:-1];
    end
    
    for j = Jindex
        for n = 1:nNeurons
            if(p1(n)>thres)
                flags(n)=1;
                CondsChange{n}=[];
                continue
            end
            
                CondsChange{n}=j;   
        end
    
        ArtifactNew = [];

        for e= 1:E
            [e CondsChange{e}]
            Gibbs       = SampleConditionalArtifact(Gibbs,e,CondsChange{e},setdiff([1:J],CondsChange{e}));
            ArtifactNew = [ArtifactNew Gibbs.variables.ArtifactE{e}];
   
   
        end
        Artifact                 = ArtifactNew;
        Gibbs.variables.Artifact = ArtifactNew;
        Gibbs                    = GibbsSamplerSpikesArtifact(Gibbs);

    end
    nIter = nIter+1;
end
   

