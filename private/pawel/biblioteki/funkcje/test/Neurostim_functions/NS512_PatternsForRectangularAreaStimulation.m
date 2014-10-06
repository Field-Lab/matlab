function Electrodes=NS512_PatternsForRectangularAreaStimulation(RowsIndexes,ColumnsIndexes);

Indexes=NS512_RowAndColumnForElectrode(500,[1:512]);

electrodes=[];
for i=1:length(RowsIndexes)
    RowIndex=RowsIndexes(i);

    for j=1:length(ColumnsIndexes)
        Column=ColumnsIndexes(j);
        electrode=find(Indexes(:,2)==RowIndex & Indexes(:,1)==Column);
        electrodes=[electrodes electrode];
    end    
    Electrodes=electrodes;
end