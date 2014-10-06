electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
figure(301)
clf;

MarkedElectrodes=[17 19];
MarkedElectrodes=PrimaryElectrodes
MovieNumbers=[14 15 17 9 12 23 12 15 15 16 11 18 10 6 18 13 9 17];

%MarkedElectrodes=ElectrodesCombinedOrder
for i=1:512    
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
    %plot(X(i),Y(i),'bd');        
end
plot(X,Y,'bd');
axis([-1000 1000 -500 500])

for i=1:length(MarkedElectrodes)
    el=MarkedElectrodes(i)
    X(i)=electrodeMap.getXPosition(el);
    Y(i)=electrodeMap.getYPosition(el);
    h=text(X(i),Y(i),num2str(el));
    %j=find(MarkedElectrodes==i);
    set(h,'FontSize',12);    
        set(h,'Color','r');
        set(h,'FontSize',18);      
    %pause(1);
    %refresh;
end
axis([-1000 1000 -500 500])