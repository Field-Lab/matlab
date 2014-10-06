ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

%cd /Volumes/Lee/Data/Lauren/2008-08-27-4/;
%cd /Volumes/Lee/Data/Lauren/2008-11-10-3/;
%cd /Volumes/Lee/Data/Lauren/2008-08-26-0/;
%cd /tmp/lhruby/Data/2008-08-26-0/
cd /Volumes/Palace/Data/Lauren/2008-08-26-0/


FileName='003';
WritePath='/Volumes/Palace/Analysis/Lauren/2008-08-26-0/data003';
%WritePath = '/snle/lab/Experiments/Array/Analysis/2008-08-26-0/data008';
%WritePathFigs='/Volumes/Lee/Analysis/Lauren/2008-08-27-4/data005Figs';
WritePathFigs = '';

FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-50 50],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants);

% 
% 
% FileName='021';
% WritePath='/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data021';
% WritePathFigs='/Volumes/Lee/Analysis/Lauren/2008-11-10-3/data021Figs';
% 
% 
% FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-50 50],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
% [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants);
% 
% 
% 
% 
% 
% cd /Volumes/Lee/Data/Lauren/2008-11-12-3/;
% FileName='001';
% WritePath='/Volumes/Lee/Analysis/Lauren/2008-11-12-3/data001';
% WritePathFigs='/Volumes/Lee/Analysis/Lauren/2008-11-12-3/data001Figs';
% 
% 
% FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-50 50],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
% [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants);
% 
% 
% 
% FileName='005';
% WritePath='/Volumes/Lee/Analysis/Lauren/2008-11-12-3/data005';
% WritePathFigs='/Volumes/Lee/Analysis/Lauren/2008-11-12-3/data005Figs';
% 
% FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-50 50],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
% [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,WritePathFigs,FigureProperties,NS_GlobalConstants);
% 
