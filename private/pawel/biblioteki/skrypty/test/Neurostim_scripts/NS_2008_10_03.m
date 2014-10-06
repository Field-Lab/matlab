NS_GlobalConstants=NS_GenerateGlobalConstants(61);
SPfilename='D:\2008-08-26-0\pattern008';
MovieNumber=28;
PatternNumber=97;
Channels=[1:7 64];
ArrayID=1;
[Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,MovieNumber,NS_GlobalConstants);
Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,PatternNumber);
[patternTimes,traces] = plotStimulusTraces(Pattern, Status, Channels, [0 50], NS_GlobalConstants);
FigureProperties=struct('FigureNumber',14,'Subplot',[2 3 3],'TimeRange',[0 50],'AmplitudeRange',[-2 2],'FontSize',16,'Colors',['k' 'r' 'b' 'm' 'g' 'c' 'y'],'LineWidth',1,'YLabel','current [\muA]');
figHandle = NS_PlotClustersOfSignaturesOnArrayLayout(traces,Channels,0,ArrayID,FigureProperties,NS_GlobalConstants,patternTimes);
Patterns(6309).data
