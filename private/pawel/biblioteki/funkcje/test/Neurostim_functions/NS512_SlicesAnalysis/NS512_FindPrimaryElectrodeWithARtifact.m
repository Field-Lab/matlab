function [PrimaryElectrode,Spikes,EI]=NS512_FindPrimaryElectrodeWithARtifact(RawDataPath,SpikeTimes,ElectrodesToCancel,Latencies);
%Finds primary electrode for specific neuron in the stimulated data. 
%Input:
%RawDataPath - path to raw data, typically after subtraction of the
%TTX-based adrtifact, but this is not mandatory
%SpikeTimes - spike times
%Patterns - the same length as SpikeTimes. Defines the latest pattern
%stimulated before the time of given spike. 
%Latencies - the same size as Spike Times, defines latency related to time
%of latest stimulation pulse
N=100;
offset=10;
LBack=20;
L=40;
spikes2=zeros(N,512,L);

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(RawDataPath); 
for i=1:N
    t=SpikeTimes(offset+i);
    d0=rawFile.getData(t-LBack,L)';
    %Electrodes0=ElectrodesToCancel(i,:);
    %d0(nonzeros(Electrodes0),:)=0;
    d1=d0([1:512]+1,:);
    spikes2(i,:,:)=d1;
end
EI=mean(spikes2);
%plot(max(abs(EI),[],3))
a=max(abs(EI),[],3);
[a1 a2]=max(a);
PrimaryElectrode=a2;
Spikes=reshape(spikes2(:,PrimaryElectrode,:),N,L);