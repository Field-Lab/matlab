ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;

RawDataPath='E:\2008-08-26-0\';
RawFileNumber='008';
ClusterFilePath='E:\analysis\2008-08-26-0\data008_proba4\protected\clusters008';
Patterns=[73:78 85:90]+110;
Movie=44;
[SuccIfSucc,SuccIfFail,PrevSucc,PrevFail]=NS_FindConditionalEfficacy(RawDataPath,RawFileNumber,ClusterFilePath,Patterns,Movie,NS_GlobalConstants);
SuccIfSucc
SuccIfFail
PrevSucc
PrevFail