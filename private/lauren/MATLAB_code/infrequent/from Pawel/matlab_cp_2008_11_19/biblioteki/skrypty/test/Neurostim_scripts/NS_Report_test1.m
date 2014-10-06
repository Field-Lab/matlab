ChipAddresses=[31 30];
NumberOfChannelsPerChip=32;

cd C:/data;
filename='004';
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
unit=1;
Fs=20000;

NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

[CI,MI]=NS_Report(filename,unit,NS_GlobalConstants);