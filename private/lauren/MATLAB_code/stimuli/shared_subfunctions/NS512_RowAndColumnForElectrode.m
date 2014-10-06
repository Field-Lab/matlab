function Indexes=NS512_RowAndColumnForElectrode(ArrayID,Channels);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

for i=1:length(Channels)
    Coordinates(i,1)=electrodeMap.getXPosition(Channels(i));
    Coordinates(i,2)=electrodeMap.getYPosition(Channels(i));
end

Indexes=Coordinates;
Indexes(:,1)=round((Coordinates(:,1)+980)./60);
Indexes(:,2)=round((Coordinates(:,2)+490)./60);