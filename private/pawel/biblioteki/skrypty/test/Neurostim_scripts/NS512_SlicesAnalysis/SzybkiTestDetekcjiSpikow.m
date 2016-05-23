DataPath='D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\data001preproc'
ArtifactDataPath=DataPath;
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,0,62,3,1500,0);

NeuronID=153;
fid = fopen('D:\Home\Pawel\analysis\slices\2013\2013-12-12-3-PH\NeuronsPatternsPrimaryElectrodes.bin','r');
Data0=fread(fid,'int32');
fclose(fid);
Data=reshape(Data0,3,length(Data0)/3);

Events=NS512_DetectSpikes(DataTraces,50,10);
[a,b]=find(Events>0);

figure(6)
n=hist(b,[1:512])

length(find(b==11))

figure(4)
plot(reshape(DataTraces(:,11,:),50,400)')

