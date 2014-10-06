function [MinValues,MaxValues,EI]=NS512_FindEIamplitudes(RawDataPath,VisionOutputPath,DataID,NeuronID);

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(RawDataPath); 

paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile([VisionOutputPath '\data' DataID '.params']);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile([VisionOutputPath '\data' DataID '.neurons'])
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

NumberOfSpikesToRead=min(length(spikeTimes),100);
L=140;
spikes=zeros(NumberOfSpikesToRead,512,L);
for i=1:NumberOfSpikesToRead
    t=spikeTimes(i);
    d0=rawFile.getData(t-20,L)';
    spikes(i,:,:)=d0(2:513,:);
end
MinValues=min(mean(spikes),[],3);
MaxValues=max(mean(spikes),[],3);
EI=mean(spikes);
%Amplitudes=MaxValues-MinValues;