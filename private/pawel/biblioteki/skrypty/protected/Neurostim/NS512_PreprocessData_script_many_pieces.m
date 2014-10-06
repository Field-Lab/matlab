ChipAddresses=[24:31];
ChipAddresses=[31:-1:24];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);

cd I:\data\2010-09-14-0

FileName='009';
WritePath = 'G:\analysis\slices\2010-09-14-0\Data_proc\data009';
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePath,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:202],[1:512],[1:512],0);

FileName='002';
WritePath = 'G:\analysis\slices\2010-09-14-0\Data_proc\data002';
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePath,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:202],[1:512],[1:512],0);



break

cd J:\2011-06-29-1

FileName='001';
WritePath = 'G:\analysis\slices\2011-06-29-1\Data_proc\data001';
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePath,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:202],[1:512],[1:512],0);

FileName='003';
WritePath = 'G:\analysis\slices\2011-06-29-1\Data_proc\data003';
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePath,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:202],[1:512],[1:512],0);


cd I:\2011-07-01-1

FileName='001';
WritePath = 'G:\analysis\slices\2011-07-01-1\Data_proc\data001';
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePath,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:240],[1:512],[1:512],0);

FileName='003';
WritePath = 'G:\analysis\slices\2011-07-01-1\Data_proc\data003';
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePath,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:240],[1:512],[1:512],0);
