function [DataTraces]=NS_ShowPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,Patterns,PatternsIndexes,FigureProperties,NS_GlobalConstants);
%This function shows the data on the array layout. 

[DataTraces,ArtifactDataTraces,Channels]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,Patterns,PatternsIndexes,FigureProperties,NS_GlobalConstants);
y=NS_AddClusterOfSignaturesToPlotOnArrayLayout(DataTraces,Channels,'b',1,FigureProperties,NS_GlobalConstants);