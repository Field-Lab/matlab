ChipAddresses=[24:31];
ChipAddresses=[31:-1:24];
NumberOfChannelsPerChip=64;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=500;

cd E:\2012-09-18-1;
FileName='001';
WritePath='C:\home\data\2012-09-18-1';
WritePathFigs=WritePath;

AmplitudeRange=[-100 100];
FigureProperties=struct('FigureNumber',1,'Subplot',[2 3 3],'TimeRange',[0 40],'AmplitudeRange',[-100 100],'FontSize',20,'Colors',['g' 'r' 'b' 'm' 'k'],'LineWidth',1);
%for i=4
    %FileName=['00' num2str(i-1)]
    [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,0,[1:14],[1:512],[1:512]);
%end