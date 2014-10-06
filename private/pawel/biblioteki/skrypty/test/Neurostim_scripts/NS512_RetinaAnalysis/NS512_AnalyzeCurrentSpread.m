electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
MovieNumber=5;
DataPath='C:\pawel\nauka\analiza\retina\2012-09-27-4\ScanCurrentSpread';

PPs=zeros(512,512);
%Distances=zeros(512,512);

for i=1:512
    Coordinates(1,i)=electrodeMap.getXPosition(i);
    Coordinates(2,i)=electrodeMap.getYPosition(i);
end

%Distances=sqrt(Coordinates(1,:).^2+Coordinates(2,:).^2);
figure(MovieNumber+10)
clf
Strange=[];
for i=1:3
    Distances(i,:)=sqrt((Coordinates(1,:)-Coordinates(1,i)).^2+(Coordinates(2,:)-Coordinates(2,i)).^2);
    PatternNumber=i
    [DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);
    Amplitude=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],PatternNumber,MovieNumber,NS_GlobalConstants);

    DataTraces=reshape(mean(DataTraces0),512,140);
    PP=max(DataTraces,[],2)-min(DataTraces,[],2);
    PPs(i,:)=PP;
    a=find(Distances(i,:)>536 & Distances(i,:)<537);
    b=find(PP>10);
    c=intersect(a,b)
    
    %subplot(5,6,i)
    loglog(Distances(i,c)',PP(c),'bd');
    loglog(Distances(i,:)',PP,'bd');
    axis([30 3000 1 1000])
    grid on
    hold on
    NeighborElectrodes=electrodeMap.getAdjacentsTo(PatternNumber,1);
    Amplitudes1(i)=mean(PP(electrodeMap.getAdjacentsTo(PatternNumber,1)));
    Amplitudes2(i)=mean(PP(electrodeMap.getAdjacentsTo(PatternNumber,2)));
    Amplitudes3(i)=mean(PP(electrodeMap.getAdjacentsTo(PatternNumber,3)));
    Amplitudes4(i)=mean(PP(electrodeMap.getAdjacentsTo(PatternNumber,4)));
end
