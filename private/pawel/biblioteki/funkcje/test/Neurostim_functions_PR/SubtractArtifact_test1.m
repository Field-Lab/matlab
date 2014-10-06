DataPath= 'D:\Home\Rydygier\Neuro\files';
PatternNumber=1;
MovieNumber=89;

[DataTraces0,ArtifactDataTraces,DataChannels]=NS_ReadPreprocessedData(DataPath,DataPath,0,PatternNumber,MovieNumber,0,0);  
[ Artifact ] = ArtifactEstimation_PR(DataTraces0,10);
h=Traces(:,121,:);
plot(reshape(h,101,80)')