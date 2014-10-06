function y=NS512_PlotGraphWithVelocities(Pairs,ArrayID,Velocities);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

clf;

for i=1:512
    x=electrodeMap.getXPosition(i);
    y=electrodeMap.getYPosition(i);
    h=plot(x,y,'ko');
    set(h,'MarkerSize',2);
    hold on;
end
    
SP=size(Pairs);
for i=1:SP(1)
    x(1)=electrodeMap.getXPosition(Pairs(i,1));
    y(1)=electrodeMap.getYPosition(Pairs(i,1));
    x(2)=electrodeMap.getXPosition(Pairs(i,2));
    y(2)=electrodeMap.getYPosition(Pairs(i,2));
               
    h=plot(x,y);
    set(h,'LineWidth',ceil(Velocities(i)));    
    h=plot(x(1),y(1),'bo');
    set(h,'MarkerSize',16);
    h=text(x(1)-28,y(1)-20,num2str(Pairs(i,1)));
    set(h,'FontSize',18);
    
    h=plot(x(2),y(2),'rd');
    set(h,'MarkerSize',16);
    h=text(x(2)-28,y(2)-20,num2str(Pairs(i,2)));
    set(h,'FontSize',16);
    
    h=gca;
    set(h,'FontSize',16);
    
    h=xlabel('microns');
    h=ylabel('microns');
        
    hold on;    
end 