NS_GlobalConstants=NS_GenerateGlobalConstants(512);
FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 140],'AmplitudeRange',[-40 20],'FontSize',11,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
NoisyChannels=[498:501 507 508];

Events=[59 14 44 49 41 32 64 7 53 64 23 52 64 35 26 64 37 26 64 51 39];
ClusterFile='F:\analysis\2009-11-20-0\data003e\ClusterFile_003';
ClusterIndexes=NS_ReadClusterFileAll(ClusterFile);
CI=ClusterIndexes(:,:,1:100);

a=sum(CI,3);
[m,p]=find(a>100); %ktore kombinacje "pattern number" i "movie number" w pliku Clusterfile 

DataPath='F:\analysis\2009-11-20-0\data003e';
WritePathFigs='F:\analysis\2009-11-20-0\data003e\figures';
for i=1 %:length(m)
    movie=m(i);
    pattern=p(i);
    [StChannels,Amplitudes]=NS_StimulatedChannels(DataPath,pattern,movie,[1:512],NS_GlobalConstants);    
    %StChannels=StChannels-320;    
    
    Indexes=CI(movie,pattern,:);
    
    [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,pattern,movie,100,0);
    
    PlotChannels=NS_RemoveBadChannels(Channels,[321:368]);
    AnalysisChannels0=NS_RemoveBadChannels(PlotChannels,StChannels);
    AnalysisChannels=NS_RemoveBadChannels(AnalysisChannels0,NoisyChannels);        
    
    DT=DataTraces(:,PlotChannels-320,:);    
    
    artifacts=find(Indexes==1);
    responses=find(Indexes==2);
    
    artifactEI=mean(DT(artifacts,:,:));
    responseEI=mean(DT(responses,:,:));
    spikeEI=responseEI-artifactEI;
        
    SsEI=size(spikeEI);
    sei=reshape(spikeEI,SsEI(2),SsEI(3));
    v=min(min(sei));
    [el,samp]=find(sei==v);
    
    PrimaryElectrodes(i)=el;
    
    y=NS512_PlotClustersOfSignaturesOnArrayLayout(spikeEI,[369:512],[1],500,FigureProperties,NS_GlobalConstants);
    FullName=[WritePathFigs '\' 'n' num2str(i) 'p' num2str(pattern) '_m' num2str(movie)];
    hj=gcf;
    %print(hj, '-dtiff', FullName);
end      