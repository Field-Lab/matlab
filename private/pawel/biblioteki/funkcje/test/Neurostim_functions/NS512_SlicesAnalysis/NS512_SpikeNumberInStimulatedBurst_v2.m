function Indexes=NS512_SpikeNumberInStimulatedBurst_v2(dane);

SD=size(dane)
Indexes=zeros(1,SD(2));

for i=1:SD(2)   
    i
    SpikesInBurst=find(dane(1,:)==dane(1,i) & dane(2,:)==dane(2,i) & dane(3,:)==dane(3,i) & dane(5,:)==dane(5,i)) %all spikes in the same burst aas the given spike - sometimes just one spike
    SpikesTimes=dane(4,SpikesInBurst)
    [T,SpikeOrder]=sort(SpikesTimes);    
    Indexes(1,SpikesInBurst(SpikeOrder))=[1:length(SpikesTimes)];
end