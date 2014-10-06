function h=NS512_ShowEIAsCircles2(EI,ArrayID,Channels,MarkedChannels,MarkedChannels2,MarkedChannels3,XRange,YRange);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

hold on;
for i=1:length(Channels)    
    S=max(abs(EI(i,:)));
    rad=max(1,round(S/50));
    
    X=electrodeMap.getXPosition(Channels(i));
    Y=electrodeMap.getYPosition(Channels(i));
    if find((MarkedChannels==Channels(i)))
        h=plot(X,Y,'ro');
        set(h,'MarkerFaceColor','r');
    else
        %if find((MarkedChannels2==Channels(i)))
        %    h=plot(X,Y,'go');
        %    set(h,'MarkerFaceColor','g');
        %else            
            h=plot(X,Y,'bo');
            set(h,'MarkerFaceColor','b');
        %end        
    end
    set(h,'MarkerSize',rad);
    
    if find((MarkedChannels2==Channels(i)))        
        h=plot(X,Y,'gd');
        set(h,'LineWidth',2);
        set(h,'MarkerSize',10);
    end
    
    if find((MarkedChannels3==Channels(i)))        
        h=plot(X,Y,'gd');
        set(h,'LineWidth',2);
        set(h,'MarkerSize',10);
    end
    %set(h,'MarkerFaceColor','b');
end
h=gca;
set(h,'XLim',XRange);
set(h,'YLim',YRange);