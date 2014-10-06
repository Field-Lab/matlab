function h=NS512_electrode_map_animation2(Electrodes);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
axis([-1000 1000 -500 500]);
hold on;
for i=Electrodes
    X=electrodeMap.getXPosition(i);
    Y=electrodeMap.getYPosition(i);
    h=plot(X,Y,'rd');
    set(h,'MarkerSize',16);
    set(h,'MarkerFaceColor','r');
    pause(0.1);
    set(h,'MarkerSize',10);
    set(h,'MarkerEdgeColor','b');
    set(h,'MarkerFaceColor','none');
end