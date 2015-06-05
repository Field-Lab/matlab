function Gibbs = Heuristic2(Gibbs)



E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
nNeurons    = Gibbs.params.nNeurons;

thres                  = 10^-4;
maxBadconds            = 4;
maxIter                = 100;
nIter                  = 1;
maxCondsSimultaneous   = 2;
numErrorConds          = 4;
maxCondsShift          = 2;


flags=zeros(nNeurons,1);

for n=1:nNeurons

    CondsChangeold{n}=[];
    cindex1(n) = 1;
    cindex2(n) = 0;
    cindex3(n) = 0;
end



nIter=1;

while(prod(flags)==0&&nIter<maxIter)


    flags=zeros(nNeurons,1);

    p1          = Gibbs.diagnostics.Logisticsp1GoF;
    condsSorted = Gibbs.diagnostics.condSortError;
    for n=1:nNeurons
    
        condsBigError(n,:)=-sort(-condsSorted(n,1:numErrorConds))
   
        if(p1(n)>thres)
            flags(n)=1;
            CondsChange{n}=[];
            CondsChangeold{n}=[];
        continue
        end
        if(cindex3(n)>=maxCondsShift)
            CondsChange{n}=[];
            CondsChangeold{n}=[];
            flags(n)=1;
        end
        
    c1 = min(condsBigError(n,cindex1(n))+cindex3(n),J);
    c2 = min(c1+cindex2(n),J);
    
    CondsChange{n} = [c1:c2];

   
    if(isequal(unique(CondsChange{n}),unique(CondsChangeold{n})))

        
        cindex2(n) = cindex2(n)+1;

        if(cindex2(n)>maxCondsSimultaneous)

            cindex2(n) = 0;
            cindex1(n) = cindex1(n) + 1;

           
        end
        
        if(cindex1(n) > numErrorConds)
            cindex1(n) = unidrnd(numErrorConds-1);          
            cindex3(n) = cindex3(n) + 1;
            cindex2(n) = 0;

        end
        
        
    else
        cindex1(n)=1;
        cindex2(n)=0;
        cindex3(n)=0;
    end
    
    


    c1 = min(condsBigError(n,cindex1(n))+cindex3(n),J);
    c2 = min(c1+cindex2(n),J);
    CondsChange{n} = [c1:c2];
    
    CondsChangeold{n} = CondsChange{n};


    listconds{n}(nIter,:)=[c1 c2 cindex1(n) cindex2(n) cindex3(n)];

    end
    ArtifactNew=[];

    for e = 1:E
    [e CondsChange{e}]
        Gibbs       = SampleConditionalArtifactT(Gibbs,e,CondsChange{e},unique([1:max(CondsChange{e})]));
        ArtifactNew = [ArtifactNew Gibbs.variables.ArtifactE{e}];
   
   
    end
    Artifact                 = ArtifactNew;
    Gibbs.variables.Artifact = ArtifactNew;
    Gibbs                    = GibbsSamplerSpikesArtifact(Gibbs);
    nIter = nIter+1;
end
Gibbs.params.listconds =listconds;




