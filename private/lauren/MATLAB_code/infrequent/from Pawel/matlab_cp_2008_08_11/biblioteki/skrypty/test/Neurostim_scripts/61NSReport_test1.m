ChipAddresses=[31 30];
NumberOfChannelsPerChip=32;

cd C:/praca;
filename='movie003';
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
unit=1;

output=61NSReportStimulatedChannels(filename,channel,time_range,ChipAddresses,NumberOfChannelsPerChip,CurrentRanges,unit);