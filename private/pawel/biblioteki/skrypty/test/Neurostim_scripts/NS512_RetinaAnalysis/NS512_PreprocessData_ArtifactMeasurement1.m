ChipAddresses=[31:-1:24];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);


cd E:\2012-09-27-4;
FileName1='012';
FileName2='013';
FileName3='014';
FileName4='015';
WritePath='F:\analiza\retina\2012-09-27-4\files\ScanCurrentSpread';
WritePathFigs=WritePath;
for i=1:16
    PatternsToRead=[1:32]+(i-1)*32;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName1,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:7],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName1,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:8],PatternsToRead,[1:512],-1);
    
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName2,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:7],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName2,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:8],PatternsToRead,[1:512],-1);

    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName3,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:7],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName3,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:8],PatternsToRead,[1:512],-1);

    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName4,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:7],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName4,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:8],PatternsToRead,[1:512],-1);
end
break;
FileName='013';
for i=1:16
    PatternsToRead=[1:32]+(i-1)*32;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:7],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:8],PatternsToRead,[1:512],-1);
end

FileName='014';
for i=1:16
    PatternsToRead=[1:32]+(i-1)*32;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:7],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:8],PatternsToRead,[1:512],-1);
end

FileName='015';
for i=1:16
    PatternsToRead=[1:32]+(i-1)*32;
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:2:7],PatternsToRead,[1:512],0);
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[2:2:8],PatternsToRead,[1:512],-1);
end