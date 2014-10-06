function [AvgPlots, AllPlots] = getTTX_PlotData( MovieNumber,PatternNumber,ChannelRead,DataPath)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%ClusterPath='C:\Users\Anna\Desktop\KSI¥¯KI\PRAKTYKI\analiza\data003\ClusterFile_006';
ArtifactDataPath=DataPath;
ArtifactSubtraction=0;
TracesNumberLimit=100; %%(albo 50)
EventNumber=0;

% %%- wczytywanie wynikow klastrowania: 26 movies,   64 patterns   53 powtorzenia
              
[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,TracesNumberLimit,EventNumber);
 DataTraces=DataTraces(1:23,:,1:40);
 SDT=size(DataTraces);
 ClusterWithArtifacts=[1:SDT(1)]';               
                
 ArtifactsMatrix=size(ClusterWithArtifacts);
 Przebieg1=DataTraces(ClusterWithArtifacts,ChannelRead,:);
 Przebieg1=(reshape(Przebieg1,ArtifactsMatrix(1),40)');
 Avg_Przebieg1_fixed=mean(Przebieg1');        

 AvgPlots=Avg_Przebieg1_fixed;
 AllPlots=Przebieg1;
end

