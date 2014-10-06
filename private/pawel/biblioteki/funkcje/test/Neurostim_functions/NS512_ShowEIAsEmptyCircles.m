function h=NS512_ShowEIAsEmptyCircles(EI,ArrayID,Channels,MarkedChannels,MarkedChannels2,XRange,YRange,Style);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

hold on;
for i=1:length(Channels)    
    S=max(EI(i,:))-min(EI(i,:));
    %MaxAbs=max(S)
    %S=S/max(S)*100;
    rad=max(1,round(S/6));
    
    X=electrodeMap.getXPosition(Channels(i));
    Y=electrodeMap.getYPosition(Channels(i));
    
    h=plot(X,Y,Style);
    set(h,'LineWidth',2);
    if find((MarkedChannels==Channels(i)))
        %h1=text(X-50,Y,num2str(Channels(i)));
        %set(h1,'Color','k');
        %set(h1,'FontSize',16);
        set(h,'MarkerFaceColor','b');
    else
        %if find((MarkedChannels2==Channels(i)))
        %    h=plot(X,Y,'go');
        %    set(h,'MarkerFaceColor','g');
        %else            
            %h=plot(X,Y,Style);
            %set(h,'LineWidth',2);
            %set(h,'MarkerFaceColor','b');
        %end        
    end
    set(h,'MarkerSize',rad);       
    
    if find((MarkedChannels2==Channels(i)))        
        h=text(X-50,Y,num2str(Channels(i)));
        dfghdhfghfgj=Style(1)
        set(h,'Color',Style(1));
        set(h,'FontSize',16);
        %h=plot(X,Y,'gd');
        %set(h,'LineWidth',2);
        %set(h,'MarkerSize',30);
    end
    %set(h,'MarkerFaceColor','b');
end
h=gca;
set(h,'XLim',XRange);
set(h,'YLim',YRange);

    
