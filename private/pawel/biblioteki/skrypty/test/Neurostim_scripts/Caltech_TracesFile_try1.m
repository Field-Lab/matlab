DataPath='C:\pawel\nauka\Caltech\wyjazd2013maj\analysis\2013-05-18'
StimulatingElectrodeNumber=2;
MovieNumber=1;

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,StimulatingElectrodeNumber,MovieNumber,0,0);

ElectrodeToLookAt=10;
Traces=reshape(DataTraces(:,ElectrodeToLookAt,:),12,2000);
figure(1)
plot(Traces')