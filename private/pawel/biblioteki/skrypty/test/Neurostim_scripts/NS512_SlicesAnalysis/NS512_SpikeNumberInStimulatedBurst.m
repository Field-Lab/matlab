function Indexes=NS512_SpikeNumberInStimulatedBurst(dane);

SD=size(dane)
Indexes=zeros(1,SD(2));

for i=1:SD(2)   
    SpikesInBurst=find(dane(1,:)==dane(1,i) & dane(3,:)==dane(3,i) & dane(5,:)==dane(5,i))
    SpikesTimes=dane(4,SpikesInBurst);
    [T,SpikeOrder]=sort(SpikesTimes);
    size(SpikeOrder);
    Indexes(1,SpikesInBurst)=SpikeOrder;
end