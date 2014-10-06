function h=NS512_ShowEIFrameAsCircles2(EI,ArrayID,Channels,MarkedChannels,Colors1,MarkedChannels2,Colors2,XRange,YRange);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

hold on;
for i=1:length(Channels)    
    S=max(EI(i,:));    
    %rad=max(1,round(S/6));
    rad=abs(EI(i,1))/4+0.1;  
  
    X=electrodeMap.getXPosition(Channels(i));
    Y=electrodeMap.getYPosition(Channels(i));
    h=plot(X,Y,'bo');
    MarkedChannelsIndex=find(MarkedChannels==Channels(i));
    if MarkedChannelsIndex        
        Color1Index=mod(MarkedChannelsIndex-1,length(Colors1))+1;
        set(h,'MarkerEdgeColor',Colors1{Color1Index});
        set(h,'MarkerFaceColor',Colors1{Color1Index});
    else                   
            %h=plot(X,Y,'bo');
            set(h,'MarkerFaceColor','b');
        %end        
    end
    set(h,'MarkerSize',rad);   
    
    MarkedChannelsIndex=find(MarkedChannels2==Channels(i));
    if MarkedChannelsIndex       
        h=plot(X,Y,'gd');
        Color2Index=mod(MarkedChannelsIndex-1,length(Colors2))+1;       
        set(h,'MarkerEdgeColor',Colors2{Color2Index});
        set(h,'LineWidth',2);
        set(h,'MarkerSize',16);
    end    
end
h=gca;
set(h,'XLim',XRange);
set(h,'YLim',YRange);

    
