function [ AvgPlots, AllPlots ] = getPlotData( MovieNumber,ClusterPath,PatternNumber,ChannelRead )
%Ta funkcja wyznacza wyznacza przebieg artefaktu dla wybranego Movie.%  

DataPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data003';
%ClusterPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data003\ClusterFile_006';
ArtifactDataPath=DataPath;
ArtifactSubtraction=0;
TracesNumberLimit=100; %%(albo 50)
EventNumber=0;

        ClusterIndexes=NS_ReadClusterFileAll(ClusterPath); %%- wczytywanie wynikow klastrowania:  26 movies,   64 patterns   53 powtorzenia
        AllClusters=ClusterIndexes(MovieNumber,PatternNumber,:);
        ClusterWithArtifacts=find(AllClusters==1);

        [DataTraces,~,~]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);

        ArtifactsMatrix=size(ClusterWithArtifacts);
        Przebieg1=DataTraces(ClusterWithArtifacts,ChannelRead,:);
        Przebieg1=(reshape(Przebieg1,ArtifactsMatrix(1),40)');
        Avg_Przebieg1=mean(Przebieg1');
        
        AvgPlots=Avg_Przebieg1;
        AllPlots=Przebieg1;
        %plot(reshape(Przebieg1,ArtifactsMatrix(1),40)')

end

