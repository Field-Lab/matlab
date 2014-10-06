function NS_GlobalConstants=NS_GenerateGlobalConstants(ElectrodeNumber);

if ElectrodeNumber==61
    ChipAddresses=[30 31];
    NumberOfChannelsPerChip=32;
else
    ChipAddresses=[31:-1:24];
    NumberOfChannelsPerChip=64;
end
        
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);