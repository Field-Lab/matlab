ChipAddresses=[24:31];
ChipAddresses=[31:-1:24];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);

tic
cd G:\2012-09-27-4
FileName='005';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan_new';
WritePathFigs=WritePath;
for i=1:16
    PatternsToRead=[1:32]+(i-1)*32;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:63],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:64],PatternsToRead,[1:512],-1);
end
toc

cd G:\2012-09-27-4
FileName='006';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan_new';
WritePathFigs=WritePath;
for i=1:16
    PatternsToRead=[1:32]+(i-1)*32;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:63],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:64],PatternsToRead,[1:512],-1);
end
toc

cd I:\2012-09-27-4
FileName='003';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan_new';
WritePathFigs=WritePath;
for i=1:16
    PatternsToRead=[1:32]+(i-1)*32;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:63],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:64],PatternsToRead,[1:512],-1);
end
toc

cd I:\2012-09-27-4
FileName='004';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan_new';
WritePathFigs=WritePath;
for i=1:16
    PatternsToRead=[1:32]+(i-1)*32;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:63],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:64],PatternsToRead,[1:512],-1);
end
toc




break
tic
cd G:\2012-09-27-4
FileName='006';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan_new';
WritePathFigs=WritePath;
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:63],[1:512],[1:512],0);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:64],[1:512],[1:512],-1);
toc
break;

cd I:\2012-09-27-4
FileName='009';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\data009';
WritePathFigs=WritePath;

AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);


tic
cd I:\2012-09-27-4
FileName='010';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\data010';
WritePathFigs=WritePath;
M1=32;
M2=16;
for i=1:M1
    Patterns=(i-1)*M2+1:i*M2;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:129],Patterns,[1:512],0);
end
toc

tic
cd I:\2012-09-27-4
FileName='011';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\data011';
WritePathFigs=WritePath;
M1=32;
M2=16;
for i=1:M1
    Patterns=(i-1)*M2+1:i*M2;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:129],Patterns,[1:512],0);
end
toc

tic
cd G:\2012-09-27-4
FileName='005';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan';
WritePathFigs=WritePath;
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:63],[1:512],[1:512],0);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:64],[1:512],[1:512],-1);
toc

tic
cd G:\2012-09-27-4
FileName='006';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\scan';
WritePathFigs=WritePath;
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:63],[1:512],[1:512],0);
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:64],[1:512],[1:512],-1);
toc

tic
cd I:\2012-09-27-4
FileName='009';
WritePath = 'D:\Home\Pawel\analysis\2012_09\2012-09-27-4\data009';
WritePathFigs=WritePath;
[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:367],[1:512],[1:512],0);
%[PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:64],[1:512],[1:512],-1);
toc
