FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-1 1],'FontSize',20,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1);

ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

%cd /Data.noindex/Lauren/2008-08-26-0/
cd D:\2008-08-26-0;
filename_movie='movie008';
l=length(filename_movie);
SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)]

patternNumberToPlot = 305;

Movie=88;
channels=[1:64];

[Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,Movie,NS_GlobalConstants);
Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,patternNumberToPlot);
[patternTimes,traces] = plotStimulusTraces(Pattern, Status, channels, [0 50], NS_GlobalConstants);

electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
%channels=electrodeMap.getAdjacentsTo(3,1)';


figHandle = NS_PlotClustersOfSignaturesOnArrayLayout(traces,channels,0,ArrayID,FigureProperties,NS_GlobalConstants,patternTimes);
%figHandle = plotStimulusTraces(filename, number_of_PD_chunk, patternNumberToPlot, channels, ArrayID, FigureProperties, NS_GlobalConstants);

