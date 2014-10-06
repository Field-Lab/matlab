electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
figure(301)
clf;

MarkedElectrodes=[1 18 26 36 37 38 41 42 43 50 51 52 53 54 56 59];
MarkedElectrodes=[16 18 27 28 37 45 51 54 60 61];
MovieNumbers=[14 15 17 9 12 23 12 15 15 16 11 18 10 6 18 13 9 17];

%MarkedElectrodes=ElectrodesCombinedOrder
%MarkedElectrodes=PatternsToUse
for i=1:512
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
    h=text(X(i),Y(i),num2str(i));
    j=find(MarkedElectrodes==i);
    set(h,'FontSize',14);
    if j
        set(h,'Color','r');
        set(h,'FontSize',22);
    end      
end
axis([-1000 1000 -500 500]);
break
figure(302);
clf;
axis([-1000 1000 -500 500]);
hold on;
for i=MarkedElectrodes
    X=electrodeMap.getXPosition(i);
    Y=electrodeMap.getYPosition(i);
    h=plot(X,Y,'rd');
    set(h,'MarkerSize',16);
    set(h,'MarkerFaceColor','r');
    pause(1);
    set(h,'MarkerSize',10);
    set(h,'MarkerEdgeColor','b');
    set(h,'MarkerFaceColor','none');
end