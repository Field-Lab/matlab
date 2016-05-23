% jak odczytac dane:
DataPath='G:\analiza\retina\2012-09-27-4\files\scan_new';
PatternNumber=317
MovieNumber=61
[DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);
s=reshape(mean(DataTraces0),512,140);
plot(s'+370)
axis([0 140 -300 300])

% jak odczytqc wspolrzedne dla konkretnej elektrody
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500);
ElectrodeNumber=400
X=electrodeMap.getXPosition(ElectrodeNumber);
Y=electrodeMap.getYPosition(ElectrodeNumber);

% jak znalexc elektrody sasiednie:
Radius=2;
Neighbors=electrodeMap.getAdjacentsTo(ElectrodeNumber,Radius)';

% 