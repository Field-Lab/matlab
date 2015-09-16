function Conds = CondsNeuron2Electrode(input,Conds)

prefElectrodes = input.neuronInfo.prefElectrodes;
nNeurons       = input.neuronInfo.nNeurons;
E              = input.tracesInfo.E;

for e = 1:E
    CondsE{e}=[];
end

for n = 1:nNeurons
    
    for e=prefElectrodes{n}

        CondsE{e}=[CondsE{e} Conds{n}];
   end
end

for e=1:E
    CondsE{e}=setdiff(sort(CondsE{e}),1);
end
Conds = CondsE;
            
    