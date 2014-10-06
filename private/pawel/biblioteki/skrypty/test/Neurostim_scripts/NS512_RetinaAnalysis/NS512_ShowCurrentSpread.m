PatternNumber=3;
MovieNumber=3;
DataPath='C:\pawel\nauka\analiza\retina\2012-09-27-4\ScanCurrentSpread';

[DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);
Amplitude=NS_AmplitudesForPattern_512_1el(DataPath,[1:512],PatternNumber,MovieNumber,NS_GlobalConstants)

DataTraces=reshape(mean(DataTraces0),512,140);
PP=max(DataTraces,[],2)-min(DataTraces,[],2);
figure(1)
clf
%PP(351)=50;
h=NS512_ShowEIFrameAsCircles(PP*2,500,[1:512],PatternNumber,[],[],[-1005 1005],[-505 505]);