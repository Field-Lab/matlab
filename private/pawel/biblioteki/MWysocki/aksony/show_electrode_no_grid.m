electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
figure(21)
for i = 1:512
    x=electrodeMap.getXPosition(i);
    y=electrodeMap.getYPosition(i);
    text(x,y,num2str(i));
    hold on
    axis([-1000 1000 -500 500])
end