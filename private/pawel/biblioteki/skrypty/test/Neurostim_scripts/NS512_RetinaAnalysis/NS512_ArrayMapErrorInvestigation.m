DataPath='G:\D\Home\Pawel\analysis\retina\2012sept\2012-09-27-4\scan_proc';
Pattern=120;
Movie=49;

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Pattern,Movie,0,0);

channels=[208 209 213 204]
for i=1:3
    subplot(3,1,i);
    s=reshape(DataTraces(1:50,Pattern-2+i,:),50,60)
    plot(s');
end