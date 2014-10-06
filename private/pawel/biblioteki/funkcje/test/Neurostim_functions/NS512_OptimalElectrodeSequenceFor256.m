function ElectrodeSequence=NS512_OptimalElectrodeSequence();
OrderFor512=NS512_OptimalElectrodeOrder();
Indexes=NS512_RowAndColumnForElectrode(500,[1:512]);

ElectrodeSequence=[];
for i=1:512
    x=OrderFor512(i,1)
    y=OrderFor512(i,2)
    
    electrode=find(Indexes(:,1)==x & Indexes(:,2)==y)
    ElectrodeSequence=[ElectrodeSequence electrode];
end
