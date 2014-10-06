function h=NS512_ShowEIAsCirclesColorCoded(EI,ArrayID,Channels,MarkedChannels,XRange,YRange);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

hold on;
for i=1:length(Channels)   
    maxes(i)=max(abs(EI(i,:)));
    mins(i)=min(abs(EI(i,:)));
end
max_global=max(maxes)
min_global=min(mins)

%MapaKolorow=[1 0 0
%    1 1 0
%    0 1 0
%    0 1 1
%    0 0 1];
N=20
MapaKolorow=jet(N);

for i=1:length(Channels)    
    X=electrodeMap.getXPosition(Channels(i));
    Y=electrodeMap.getYPosition(Channels(i));
    if find((MarkedChannels==Channels(i)))
        h=plot(X,Y,'ro');
        set(h,'MarkerFaceColor','r');
    else
        h=plot(X,Y,'bo');
        set(h,'MarkerFaceColor','b');
    end
    set(h,'MarkerSize',25);
    ColorIndex=round((maxes(i)-min_global)/(max_global-min_global)*(N-1));
    set(h,'MarkerFaceColor',MapaKolorow(ColorIndex+1,:));
end
h=gca;
set(h,'XLim',XRange);
set(h,'YLim',YRange);

    
