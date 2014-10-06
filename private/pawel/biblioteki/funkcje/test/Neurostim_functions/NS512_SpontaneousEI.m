function EI=NS512_SpontaneousEI(RawDataPath,NeuronFilePath,NeuronID);

N=200; %number of spikes
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(RawDataPath); 
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(NeuronFilePath);
spikeTimes = neuronFile.getSpikeTimes(NeuronID)';
L=80;

for i=1:min(N,length(spikeTimes))
    t=spikeTimes(i);
    d0=rawFile.getData(t-20,L)';
    size(d0)
    d1=d0([2:513],:);
    spikes2(i,:,:)=d1;
end

EI=mean(spikes2);