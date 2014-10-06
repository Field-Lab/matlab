ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

%cd /Volumes/Lee/Data/Lauren/2008-08-27-4/; %path to raw data (don't include data*** folder)
%cd /tmp/Data/lhruby/2010-03-05-4/
cd /Volumes/Palace/Data/Lauren/2010-03-05-4/
%cd /Volumes/Palace/Data/Lauren/2009-09-03-1/

FileName='004'; %folder within current directory
%WritePath='/Volumes/Lee/Analysis/Lauren/testPlotter';
%WritePath = '/tmp/Data/lhruby/2010-03-05-4/figs';
WritePath='/Volumes/Palace/Analysis/Lauren/2010-03-05-4/figs';
%WritePath = '/Volumes/Palace/Analysis/Lauren/2009-09-03-4/figs';

FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',...
    [0 40],'AmplitudeRange',[-800 200],'FontSize',20,'Colors',...
    ['g' 'r' 'b' 'm' 'k'],'LineWidth',1);


[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,y]=NS_OverlaidResponseFigureGenerator(FileName,WritePath,FigureProperties,NS_GlobalConstants);
%[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew_plot(FileName,WritePath,FigureProperties,NS_GlobalConstants, 100);