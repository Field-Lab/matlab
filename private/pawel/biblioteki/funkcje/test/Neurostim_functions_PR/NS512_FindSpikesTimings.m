function SpikesTimings=NS512_FindSpikesTimings(Traces,Events,Channels);
%Traces must be after artifact substraction
ST=size(Traces)
SpikesTimings=zeros(ST(1),length(Channels));
for i=1:length(Channels)
    WaveformTypes=Events(:,Channels(i));
    SpikeIndexes=find(WaveformTypes==1);    
    TracesSize = size(Traces);
    ChannelTraces=Traces(SpikeIndexes,Channels(i),:);
    ChannelTraces2D = reshape(ChannelTraces,length(SpikeIndexes),TracesSize(3));
    USI = FindUnifiedSpikes_PR(ChannelTraces2D,-25); %Wektor wskaznikow (numery probek dla ktorych nastapilo przekroczenie progu)        
    SpikesTimings(SpikeIndexes,i)=USI;
end