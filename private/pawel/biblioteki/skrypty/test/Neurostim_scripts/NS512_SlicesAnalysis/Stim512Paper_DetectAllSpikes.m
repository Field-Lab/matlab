DataPath='E:\analysis\2010-09-14-0\data002_preproc';
ArtifactDataPath=DataPath;
Pattern=339;
Movie=98;
channel=456;

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,Pattern,Movie,0,0);
spikes=NS512_DetectAllSpikes(DataTraces(:,:,:),[1:512],-40.1,-20.1)





break
figure(1)
clf
hold on

plot(reshape(DataTraces(:,channel,:),50,600)','b-');

SpikesForChannel=find(spikes(:,1)==channel);
for i=1:length(SpikesForChannel)
    SpikeIndex=SpikesForChannel(i);
    TraceIndex=spikes(SpikeIndex,2)
    SampleIndex=spikes(SpikeIndex,3)
    if SampleIndex>10 & SampleIndex<590
        plot(SampleIndex-5:SampleIndex+10,reshape(DataTraces(TraceIndex,channel,SampleIndex-5:SampleIndex+10),1,16)','r-');
    end
end
grid on