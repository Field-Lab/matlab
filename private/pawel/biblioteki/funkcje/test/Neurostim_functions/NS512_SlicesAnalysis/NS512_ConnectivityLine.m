function dx=NS512_ConnectivityLine(X,Y,Color1,Color2,NumberOfSteps,StepsToPlot);
% Jesli StepsToPlot jest 0, to wtedy wszystkie kroki s? rysowane. Jesli np,
% rysujemy tylko pierwszy krok, to wtedy nie ma ca?ych linii ??cz?cych
% elektrody z neuronami, tylko "kikuty" wystaj?ce z eletrod (pozwala
% pokaza? ile takich po?acze? wychodzi  z danej elektrody).

%if StepsToPlot==0
%    StepsToPlot=[1:NumberOfSteps];
%end

for i=1:length(X)-1
    x1=X(i);
    y1=Y(i);
    x2=X(i+1);
    y2=Y(i+1);
    
    dx=(x2-x1)/NumberOfSteps;
    dy=(y2-y1)/NumberOfSteps;
    
    for j=1:NumberOfSteps
        if find(StepsToPlot==j)
            color=Color1+(Color2-Color1)*j/(NumberOfSteps);
            x=x1+(j-1)*dx;
            y=y1+(j-1)*dy;
            h=plot([x1+(j-1)*dx x1+j*dx],[y1+(j-1)*dy y1+j*dy]);            
            set(h,'Color',color);
            set(h,'LineWidth',2);            
        end
    end
end