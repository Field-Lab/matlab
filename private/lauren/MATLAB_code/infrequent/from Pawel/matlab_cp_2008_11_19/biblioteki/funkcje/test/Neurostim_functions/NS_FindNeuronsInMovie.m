function Indexes=NS_FindNeuronsInMovie(Traces,Threshold);
%Input:
%Traces - array of size NxCxT, where:
%N - number of traces for each channel;
%C - number of channels;
%T - length of single trace.
TimeRange=[0 ]; %When exceeding thresholds will be founfd for sample N, samples from N-5 to N+20 will be taken as a spike shape;
if TimeRange(1)<0
    z=zeros(1,-TimeRange(1));
else
    z=[];
end
STraces=size(Traces)
SpikeTimes=zeros(1,STraces(1));

for i=1:STraces(1) %for each event...
    SpikeTime=0;
    for j=1:STraces(2) %for each channel...
        %a=Traces(i,j,1:STraces(3)-TimeRange(1)) %resize(Traces(i,j,1:STraces(3)-TimeRange(1)),1,STraces(3))
        %STraces(3)-TimeRange(2)
        %size(a)
        %b=reshape(a,1,30)
        signal=[z reshape(Traces(i,j,1:STraces(3)-TimeRange(1)),1,STraces(3))];
        T=find(abs(signal)>Threshold);
        if length(T)>0
            if SpikeTime==0 %if none spike was detected in previously chcecked channels for this event...
                SpikeTime=min(T);
            else
                if min(T)<SpikeTime %is there were some spikes detected, but in this channel there is a spike with smaller latency...
                    SpikeTime=min(T);
                end
            end
        end
    end %OK, we have found time for a spike, if the time is 0 then this event is artifact only
    SpikeTimes(i)=SpikeTime;
    %if SpikeTime~=0
    %    Spike=Traces(i,:,SpikeTime+TimeRange(1):SpikeTime+TimeRange(2));
    %    size(Spike)
    %    Spikes=[Spikes' Spike']';
    %end
end
Indexes=SpikeTimes;