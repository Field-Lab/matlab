ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

cd /Volumes/Creampuff/Data/Lauren/2008-08-27-4/;

FileName='001';
WritePath='/Volumes/Creampuff/Analysis/Lauren/2008-08-27-4/data001/analysis_images';
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-1200 200],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);

[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,y]=NS_PlotOverlaidRespForMovieAllPatterns(FileName,WritePath,FigureProperties,NS_GlobalConstants);