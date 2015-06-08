function CondDel = detectHighSpiking(Gibbs,input,neuronIndex,breakRange)


HighThres       = input.params.Heuristic.HighThres;
ActivationThres = input.params.Heuristic.ActivationThres;


E           = Gibbs.params.E;
J           = Gibbs.params.J;
T           = Gibbs.params.T;
I           = Gibbs.params.I;
nNeurons    = Gibbs.params.nNeurons;
spikes      = Gibbs.variables.spikes;
prefElectrodes  = input.neuronInfo.prefElectrodes;
e               =  prefElectrodes{neuronIndex}(1);
breakPoints     = [input.tracesInfo.breakPoints{e} J];


for n=1:nNeurons
    CondDel{n} = [];
end

if(breakRange==1)
    breakLim(1) = 1;
    breakLim(2) = breakPoints(1);
else
    breakLim(1)  = breakPoints(breakRange-1)+1;
    breakLim(2)  = breakPoints(breakRange);
end

condsActivation     = detectActivation(input,spikes,ActivationThres);
condsHighActivation = detectActivation(input,spikes,HighThres);




if(~isequal(condsActivation(neuronIndex),condsHighActivation(neuronIndex)))
	return
else
    
    if(condsActivation(neuronIndex)<=breakLim(2)&&condsActivation(neuronIndex)>=breakLim(1))
        CondDel{neuronIndex}=condsActivation(neuronIndex);
    end
end



