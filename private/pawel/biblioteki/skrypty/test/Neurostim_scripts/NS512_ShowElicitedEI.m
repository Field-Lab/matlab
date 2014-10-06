NS_GlobalConstants=NS_GenerateGlobalConstants(512);

FilePath='D:\analysis\cultures\2009-11-20-0\data003e\ClusterFile_003_n65';
ClusterIndex=NS_ReadClusterFile(FilePath,54,63,200);

DataPath='D:\analysis\cultures\2009-11-20-0\data003e\';
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,63,54,200,0);

artifacts=find(ClusterIndex==1);
spikes=find(ClusterIndex==2);
Artifact=mean(DataTraces(artifacts,:,:));
S=mean(DataTraces(spikes,:,:));
Spike=S-Artifact;
m=mean(Spike,3);
for i=1:192
    Spike(1,i)=Spike(1,i)-m(i);
end
T(2,:,5:104)=Spike(1,105:192,1:100)*2.6;

Chns=[444 445 451 452 453 460 461];
T1=T(:,Chns-424,:);
T1=T1./0.84;
FigureProperties=struct('FigureNumber',3,'TimeRange',[0 80],'AmplitudeRange',[-80 40],'FontSize',16,'Colors',['k' 'r' 'b' 'y'],'XTick',[0 0.5 1 1.5 2 2.5 3 3.5],'YTick',[-80:20:40],'LineWidth',2,'YLabel','signal [\muV]');
y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(T1,[Chns],[1 2],500,FigureProperties,NS_GlobalConstants,[444]);
break;
h=gca;
set(h,'XTickLabel',{'0' '' '1' '' '2' '' '3' '' '4'});