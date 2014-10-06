function ElectrodeSequence=NS512_OptimalElectrodeSequenceLeftHalf();
OrderFor256=NS512_OptimalElectrodeOrderLeftHalf
Indexes=NS512_RowAndColumnForElectrode(500,[1:512]);

ElectrodeSequence=[];
for i=1:256
    x=OrderFor256(i,1)
    y=OrderFor256(i,2)
    
    electrode=find(Indexes(:,1)==x & Indexes(:,2)==y)
    ElectrodeSequence=[ElectrodeSequence electrode];
end
