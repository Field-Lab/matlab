cd C:\pawel\nauka\analiza\retina\Sept2012
load Thresholds27_4;
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);

radius=3

figure(4)
MinThresholds=ones(1,512)*100;
NumberOfElectrodesToShow=zeros(1,512);
for i=1:512
    ElectrodesToExclude=electrodeMap.getAdjacentsTo(Patterns(i),radius)';
    ThresholdsForPattern=Thresholds(i,:);
    ThresholdsForPattern(ElectrodesToExclude)=0;
    ElectrodesToShow=find(ThresholdsForPattern>0); 
    NumberOfElectrodesToShow(i)=length(ElectrodesToShow);
    if ElectrodesToShow
        MinThresholds(i)=min(ThresholdsForPattern(ElectrodesToShow));
    end
end
plot(MinThresholds,'bd')
figure(13)
clf

ChannelsToPlot=find(MinThresholds<100)
NS512_ShowEIAsCirclesColorCoded(MinThresholds(ChannelsToPlot)',500,ChannelsToPlot,[],[-1005 1005],[-505 505]);
figure(14)
clf
NS512_ShowEIAsCircles(MinThresholds(ChannelsToPlot)'.^(-1)*2000,500,ChannelsToPlot,[],[-1005 1005],[-505 505]);
break
for i=1:512  
    X(i)=electrodeMap.getXPosition(i);
    Y(i)=electrodeMap.getYPosition(i);
end
figure(7)
YPrim=(Y+450)/60;
for i=0:15
    Electrodes=find(YPrim==i)
    [q,w]=sort(X(Electrodes))
    subplot(16,1,i+1)
    plot(MinThresholds(Electrodes(w)),'bd')
    axis([1 32 13 32])
    h=gca
    %set(h,'Visible','off');
    %ElectrodesOrderInRow=sort(
end

figure(6)
YPrim=(Y+450)/60;
for i=0:15
    Electrodes=find(YPrim==i)
    [q,w]=sort(X(Electrodes))
    subplot(16,1,i+1)
    plot(MinThresholds(Electrodes(w)),'bd')
    axis([1 32 13 32])
    h=gca
    %set(h,'Visible','off');
    %ElectrodesOrderInRow=sort(
end

break
for j=1:512
    X(j)=electrodeMap.getXPosition(j);
    Y(j)=electrodeMap.getYPosition(j);
end

for i=1:512
    clear X1;
    clear Y1
    ElectrodesToExclude=electrodeMap.getAdjacentsTo(Patterns(i),radius)';
    ThresholdsForPattern=Thresholds(i,:);
    ThresholdsForPattern(ElectrodesToExclude)=0;
    ElectrodesToShow=find(ThresholdsForPattern>0); 
        
    for j=1:length(ElectrodesToShow)
        X1(j)=electrodeMap.getXPosition(ElectrodesToShow(j));
        Y1(j)=electrodeMap.getYPosition(ElectrodesToShow(j));
    end
    figure(12)
    clf
    h1=plot(X,Y,'bo')
    set(h1,'MarkerSize',1);
    hold on
    %figure(13)
    h2=plot(X1,Y1,'ro')                    
    set(h2,'MarkerFaceColor','r')
    axis([-1000 1000 -500 500]);
    Thresholds2(i,:)=ThresholdsForPattern;
    pause(0.5)
end