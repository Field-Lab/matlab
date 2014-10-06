function EI=NS512_EI_Corrected_Timings(Traces,Channels,Type,UniSpikesIndic,WaveformTypes,NumberOfSamples)
%Type: 0 - artifact, 1 - spikes, 2 - different class of spikes etc.

Indexes=find(WaveformTypes==Type); % only 
TracesIndexed=Traces(Indexes,Channels,:);
TracesForEI=zeros(length(Indexes),length(Channels),NumberOfSamples);
for i=1:length(Indexes)

    ind=UniSpikesIndic(Indexes(i));
    if ind <= 40
    TracesForEI(i,:,:)=TracesIndexed(i,:,ind:ind+NumberOfSamples-1);    %Wyrazenie poczatkowe
    elseif ind > 40 && ind <50
        TracesForEI(i,:,:)=TracesIndexed(i,:,size(TracesIndexed,3)-NumberOfSamples+1:size(TracesIndexed,3));
     else
        TracesForEI(i,:,:)=zeros(1,size(TracesIndexed,2),NumberOfSamples);
    end   
end
% ind;
 STracesForEI=size(TracesForEI);
% SizeOfMeanTraces = size(mean(TracesForEI));
EI=reshape(mean(TracesForEI,1),STracesForEI(2),STracesForEI(3));