function Vectors=NS_ConcacenateWaveforms(Traces);
%Traces - 3-dimensional: NxCXL, where N - number of traces for each
%channel, C - number of channels, L - length of each waveform

STraces=size(Traces);
Vectors=zeros(STraces(1),STraces(2)*STraces(3));
for i=1:STraces(1)
    data0=Traces(i,:,:);
    data1=reshape(data0,STraces(2),STraces(3));    
    Vectors(i,:)=reshape(data1',1,STraces(2)*STraces(3)); 
end