function [PrimaryElectrode,Amplitudes]=NS512_FindPrimaryElectrode(RawDataPath,VisionOutputPath,DataID,NeuronID);

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(RawDataPath); 

%paramsFile=edu.ucsc.neurobiology.vision.io.ParametersFile([VisionOutputPath '\data' DataID '.params']);
neuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile([VisionOutputPath '\data' DataID '.neurons'])
spikeTimes = neuronFile.getSpikeTimes(NeuronID)'; % for given neuron, import the spikes times

NumberOfSpikesToRead=min(length(spikeTimes),400);
L=60;
spikes=zeros(NumberOfSpikesToRead,512,L);
for i=1:NumberOfSpikesToRead-2
    t=spikeTimes(i+2);
    d0=rawFile.getData(t-20,L)';
    spikes(i,:,:)=d0(2:513,:);
end
MinValues=min(abs((mean(spikes))),[],3);
MaxValues=max(abs((mean(spikes))),[],3);
Amplitudes=MaxValues-MinValues;
[MaxAmplitude,PrimaryElectrode]=max(Amplitudes);