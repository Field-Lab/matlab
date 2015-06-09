function CondsActivation = detectActivation(input,spikes,thres)

I           = input.tracesInfo.I;
J           = input.tracesInfo.J;
nNeurons    = length(spikes);
for n=1:nNeurons
    CondsActivation(n) = J+1;
    countSpikes       = nansum(spikes{n}');
    spikeProbs        = countSpikes./I;
    aux               = find(spikeProbs>=thres);
    if(~isempty(aux))
        CondsActivation(n)=aux(1);
    end
end
        
