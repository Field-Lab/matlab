clear;
% 1. Files, electrode numbers and constants
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
NS_GlobalConstants=NS_GenerateGlobalConstants(61);
%DataPath='C:\Users\pawel\analysis\2008-12-06-0\data002';

DataPath='H:\analysis\2008-12-06-0\data005';
ClusterFileName=[DataPath '\ClusterFile_005'];
ArtifactDataPath=DataPath;
PatternNumber=44;
StimElectrode=PatternNumber;
CenterChannel=StimElectrode;
CenterChannel=44; %may be different than the StimElectrode
Radius=0;
ArtifactSubtraction=0;
EventNumber=0; %all events - basically do not change
Channels=electrodeMap.getAdjacentsTo(CenterChannel,Radius)';

% 2. Data reading
StartMovie=39;

% 2a) Artifact and amplitude for current 1
MovieNumber=StartMovie;
FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
load(FullName_status);
[DataTraces1,ArtifactDataTraces,Channels2]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,100,EventNumber);
SD=size(DataTraces1);
WaveformTypes=NS_ReadClusterFile(ClusterFileName,MovieNumber,PatternNumber,SD(1));
a=find(WaveformTypes==1);
DataTraces1=DataTraces1(a,Channels,:);
Artifact1=mean(DataTraces1);
FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
load(FullName_pattern);
[patternTimes,traces] = plotStimulusTraces(Pattern, Status, Channels, [0 50], NS_GlobalConstants);
amplitude1=max(max(traces));

% 2b) Artifact and amplitude for current 2
MovieNumber=StartMovie+1;
FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
load(FullName_status);
[DataTraces1,ArtifactDataTraces,Channels2]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,100,EventNumber);
SD=size(DataTraces1);
WaveformTypes=NS_ReadClusterFile(ClusterFileName,MovieNumber,PatternNumber,SD(1));
a=find(WaveformTypes==1);
DataTraces1=DataTraces1(a,Channels,:);
Artifact2=mean(DataTraces1);
FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
load(FullName_pattern);
[patternTimes,traces] = plotStimulusTraces(Pattern, Status, Channels, [0 50], NS_GlobalConstants);
amplitude2=max(max(traces));

% 2c) The real data
MovieNumber=StartMovie+2;
FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
load(FullName_status);
[DataTraces,ArtifactDataTraces,Channels2]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,100,EventNumber);
SD=size(DataTraces1);
WaveformTypes=NS_ReadClusterFile(ClusterFileName,MovieNumber,PatternNumber,SD(1));
DataTraces=DataTraces(:,Channels,:);
FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
load(FullName_pattern);
[patternTimes,traces] = plotStimulusTraces(Pattern, Status, Channels, [0 50], NS_GlobalConstants);
amplitude3=max(max(traces));

Artifact3=Artifact2+(amplitude3-amplitude2)/(amplitude2-amplitude1)*(Artifact2-Artifact1);
SD=size(DataTraces);
for i=1:SD(1)
    DataTraces(i,:,:)=DataTraces(i,:,:)-Artifact3;
end
%break;
% * * * plotting
LowLimit=-400;
HighLimit=400;
FigureProperties=struct('FigureNumber',22,'Subplot',[2 3 3],'TimeRange',[5 35],'AmplitudeRange',[LowLimit HighLimit],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','signal [mV]');
y=NS_PlotClustersOfSignaturesOnArrayLayout(DataTraces,Channels,WaveformTypes,1,FigureProperties,NS_GlobalConstants);
%y=NS_PlotClustersOfSignaturesOnArrayLayout(Artifact1,Channels,[1],1,FigureProperties,NS_GlobalConstants);

StartMovie=6;
for i=1:10
    MovieNumber=StartMovie+i*3;
    FullName_status=[DataPath filesep 'status_m' num2str(MovieNumber)];
    load(FullName_status);
    [DataTraces1,ArtifactDataTraces,Channels2]=NS_ReadPreprocessedData(DataPath,ArtifactDataPath,ArtifactSubtraction,PatternNumber,MovieNumber,100,EventNumber);
    SD=size(DataTraces1);
    WaveformTypes=NS_ReadClusterFile(ClusterFileName,MovieNumber,PatternNumber,SD(1));
    a=find(WaveformTypes==1);
    DataTraces1=DataTraces1(a,Channels,:);
    Artifact1=mean(DataTraces1);
    FullName_pattern=[DataPath filesep 'pattern' num2str(PatternNumber) '_m' num2str(MovieNumber)];
    load(FullName_pattern);
    [patternTimes,traces] = plotStimulusTraces(Pattern, Status, Channels, [0 50], NS_GlobalConstants);
    amplitude1=max(max(traces));
end

    