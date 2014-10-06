function Patterns=NS512_PatternsForRowStimulation(RowIndexes,WhichColumns);

Indexes=NS512_RowAndColumnForElectrode(500,[1:512]);
Patterns=zeros(length(RowIndexes),length(WhichColumns));

for i=1:length(RowIndexes)
    RowIndex=RowIndexes(i);
    electrodes=[];
    for j=1:length(WhichColumns)
        Column=WhichColumns(j);
        electrode=find(Indexes(:,2)==RowIndex & Indexes(:,1)==Column);
        electrodes=[electrodes electrode];
    end    
    Patterns(i,:)=electrodes;
end