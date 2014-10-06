function Angles=NS512_ConnectivityAngle(Pairs,ArrayID);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);
SP=size(Pairs);

for p=1:SP(1)
    x1=electrodeMap.getXPosition(Pairs(p,1));
    y1=electrodeMap.getYPosition(Pairs(p,1));
    x2=electrodeMap.getXPosition(Pairs(p,2));
    y2=electrodeMap.getYPosition(Pairs(p,2));
    
    a1=x1+sqrt(-1)*y1;
    a2=x2+sqrt(-1)*y2;
    
    Angles(p)=angle(a1-a2)*360/6.28;
end
 