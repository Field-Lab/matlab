function [f,spectrum]=NS_CalculateMultipleFFTs(Traces,SamplingFrequency);

N=4096;
STraces=size(Traces);
f=[0:STraces(2)-1]/STraces(2)*SamplingFrequency;

f=[0:N-1]/N*SamplingFrequency;
spectrum=zeros(STraces(1),N);

for i=1:STraces
    s=Traces(i,:);
    spectrum(i,:)=fft(s,N);
end