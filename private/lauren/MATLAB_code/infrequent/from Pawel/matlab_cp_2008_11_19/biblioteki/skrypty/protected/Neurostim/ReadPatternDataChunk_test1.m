cd C:\praca\stimuli_files;
SPfilename='pattern000';
number_of_PD_chunk=2;
ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;

ChipAddresses=[30 31];
NumberOfChannelsPerChip=32;
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040];
Fs=20000;
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);
ArrayID=1;
tic
[patterns_out,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,number_of_PD_chunk,NS_GlobalConstants);
toc
for i=1:0
    if Status.ChannelsStatus(i).range~=2
        i
        Status.ChannelsStatus(i).range
    end
end

for i=1:0
    name=NS_PatternAmplitudes(patterns_out,PatternsIndexes,Status,i,NS_GlobalConstants)
end

%Pulse=NS_FindPulsesShapesForPattern('008',3,2,NS_GlobalConstants,patterns_
%out);