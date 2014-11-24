function select_clusters(cones)

global myMap conesSavePosition myCluster

hPlot=findobj('position',conesSavePosition);
subplot(hPlot);
hold on
myXaxis=get(gca,'XLim');
myYaxis=get(gca,'YLim');
    
allCones=cell2mat(cones');

flag=1;
cnt=1;
if isempty(myCluster)
    myCluster=cell(20,1);    
end

i=1;
while ~isempty(myCluster{i})
    i=i+1;
end

markers='*vd>+<x.........';
colors=[1 1 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0; 1 1 1; 1 0 0; 0 1 0; 0 1 1; 1 0 1; 1 1 0];

t=unique(myMap);
for j=1:length(t)
    [pos_y(j),pos_x(j)]=find(myMap==t(j),1);
end
while flag
    

    coords=round(getrect);
    x=[coords(1) coords(1)+coords(3)];
    y=[coords(2) coords(2)+coords(4)];
    
    if x(1)>myXaxis(1) && x(2)<myXaxis(2) && y(1)>myYaxis(1) && y(2)<myYaxis(2)        
        
        polyg_x=[x(1), x(1), x(2), x(2), x(1)];
        polyg_y=[y(1), y(2), y(2), y(1), y(1)];        
        
        IN = inpolygon(pos_x,pos_y,polyg_x,polyg_y);        
        
        myCluster{i}=[myCluster{i}; t(IN)];
        plot(pos_x(IN),pos_y(IN),'color',colors(i,:), 'marker', markers(i))

    else
        flag=0;
    end

    
end
