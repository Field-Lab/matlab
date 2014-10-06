function NeuronEI=NS512_LoadEIFromRawData(RawDataPath,spikeTimes,Length,Offset);
%Ta funkcja laduje dane z plikow RAW, bazujac na podanych spikeTimes.
%Generalnie, moze byc rownie dobrze uzyta do zaladowania np. EI neuronu,
%lub tez spike-triggered EI.

rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(RawDataPath); 

spikes=zeros(513,Length);
for i=1:length(spikeTimes)
    t=spikeTimes(i)-Offset;
    d0=rawFile.getData(t,Length)';
    spikes=spikes+double(d0);
end
NeuronEI=spikes(2:513,:)/length(spikeTimes);