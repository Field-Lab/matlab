NS_GlobalConstants=NS_GenerateGlobalConstants(512);
ChipAddresses=[24:31];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
%NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

ArrayID=500;

cd I:\2009-11-27-0;
FileName='001';
WritePath='E:\pawel\analysis\retina\2009-11-27-0\data001\April2010\files';

WritePathFigs=WritePath;

FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-50 50],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,0,[],[],[1:512]);

clear;
Ns512_DetectSpikes_test2;
break;

FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-50 50],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
el=[20 70 76 106 139 354 358 383 421 444 445 463];
%el=[16 25 26 36 65 81 100 105 120 144];
electrodes=[];
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(500); %define the electrode map - must be different than 1 for silicon probes
Radius=1;
for i=el
    Channels=electrodeMap.getAdjacentsTo(i,Radius)';
    electrodes=[electrodes Channels];
end
electrodes=unique(electrodes);


%break
electrodes=424:468;
WritePath='D:\analysis\2009-11-20-0\data003c';
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,0,[1:2],[1:5],electrodes);
break;
electrodes=469:512;
WritePath='D:\analysis\2009-11-20-0\data003d';
%[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants,0,[1:2],[1:5],electrodes);

PrepDataPath='D:\analysis\retina\2009-11-27-0\data001';
ClusterFileName='D:\analysis\retina\2009-11-27-0\data001\ClusterFile_001';
GoodChannels=[385:512];
Movies=[1:2:151];

Patterns=[1:64];
Responses=zeros(length(Patterns),512);
Noise=Responses;
tic
for PatternNumber=Patterns
    PatternNumber
    [a,b]=NS512_ApplyLinearArtifactModel2(PrepDataPath,PatternNumber,Movies,GoodChannels,0,0,ClusterFileName,NS_GlobalConstants);
    Responses(PatternNumber,:)=a;
    Noise(PatternNumber,:)=b;
end
toc
save Responses;
Save Noise;