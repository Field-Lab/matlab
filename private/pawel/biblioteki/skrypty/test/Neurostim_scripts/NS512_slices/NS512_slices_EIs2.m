NS_GlobalConstants=NS_GenerateGlobalConstants(512);
FigureProperties=struct('FigureNumber',12,'Subplot',[2 3 3],'TimeRange',[0 140],'AmplitudeRange',[-20 10],'FontSize',11,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
NoisyChannels=[498:501 507 508];
Electrode0=320;

Events=[59 14 44 49 41 32 64 7 53 64 23 52 64 35 26 64 37 26 64 51 39];
ClusterFile='D:\analysis\2009-11-20-0\data003e\ClusterFile_003';
ClusterIndexes=NS_ReadClusterFileAll(ClusterFile);
CI=ClusterIndexes(:,:,1:100);

a=sum(CI,3);
[m,p]=find(a>100); %ktore kombinacje "pattern number" i "movie number" w pliku Clusterfile 

DataPath='D:\analysis\2009-11-20-0\data003e';
WritePathFigs='D:\analysis\2009-11-20-0\data003e\FiguresAndMovies';
for i=1:length(m)
    movie=m(i);
    pattern=p(i);
    [StChannels,Amplitudes]=NS_StimulatedChannels(DataPath,pattern,movie,[1:512],NS_GlobalConstants);    
    %StChannels=StChannels-320;    
    
    Indexes=CI(movie,pattern,:);    
    artifacts=find(Indexes==1);
    responses=find(Indexes==2);
    
    [DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,DataPath,0,pattern,movie,100,0);
    
    %1. analiza
    AnalysisDataTraces=DataTraces;
    for j=NoisyChannels
        AnalysisDataTraces(:,j-Electrode0,:)=0;
    end
    
    for j=StChannels
        AnalysisDataTraces(:,j-Electrode0,:)=0;
    end
                                       
    artifactEI=mean(AnalysisDataTraces(artifacts,:,:));
    responseEI=mean(AnalysisDataTraces(responses,:,:));
    spikeEI=responseEI-artifactEI;        
    SsEI=size(spikeEI);
    sei=reshape(spikeEI,SsEI(2),SsEI(3));
    v=min(min(sei));
    [el,samp]=find(sei==v);    
    PrimaryElectrodes(i)=el;
    
    %DT=DataTraces(:,PlotChannels-320,:);    
    artifactEI=mean(DataTraces(artifacts,49:192,:));
    responseEI=mean(DataTraces(responses,49:192,:));
    spikeEI=responseEI-artifactEI;
    
    %y=NS512_PlotClustersOfSignaturesOnArrayLayoutWithMarks(spikeEI,[369:512],[1],500,FigureProperties,NS_GlobalConstants,StChannels);
    
    %h=gcf;
    %FullName=[WritePathFigs '\' 'n' num2str(i) 'p' num2str(pattern) '_m' num2str(movie)];    
    %hj=gcf;
    %print(hj, '-dtiff', FullName);
    %Mark the stimulation pulse
    for j=1:length(StChannels)
        k=find(Channels==StChannels(j));
        sei(k,[1 2 5 6])=max(max(abs(sei)))/3;
        sei(k,[3 4])=-max(max(abs(sei)))/3;
    end
    
    %M=NS512_SaveMovieFromSignatureWithMarks(sei(:,1:100),Channels,StChannels,500,FigureProperties,NS_GlobalConstants,StChannels);        
    %OutputMovie=NS512_ReplicateMovieFrames(M,4);   
    %mpgwrite(OutputMovie, 'jet', FullName);
end      