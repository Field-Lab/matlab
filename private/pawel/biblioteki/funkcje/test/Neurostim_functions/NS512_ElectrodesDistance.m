function Distance=NS512_ElectrodesDistance(ElectrodeFrom,ElectrodesTo,ArrayID);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

xFrom=electrodeMap.getXPosition(ElectrodeFrom);
yFrom=electrodeMap.getYPosition(ElectrodeFrom);

for i=1:length(ElectrodesTo)
    x=electrodeMap.getXPosition(ElectrodesTo(i));
    y=electrodeMap.getYPosition(ElectrodesTo(i));    
    Distance(i)=sqrt((xFrom-x)^2+(yFrom-y)^2);
end