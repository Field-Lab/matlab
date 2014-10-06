function Patterns=NS512_PatternsForColumnStimulation(ColumnIndexes,WhichRows);

Indexes=NS512_RowAndColumnForElectrode(500,[1:512]);
Patterns=zeros(length(ColumnIndexes),length(WhichRows));

for i=1:length(ColumnIndexes)
    ColumnIndex=ColumnIndexes(i);
    electrodes=[];
    for j=1:length(WhichRows)
        Row=WhichRows(j);
        electrode=find(Indexes(:,2)==Row & Indexes(:,1)==ColumnIndex);
        electrodes=[electrodes electrode];
    end        
    Patterns(i,:)=electrodes;
end


