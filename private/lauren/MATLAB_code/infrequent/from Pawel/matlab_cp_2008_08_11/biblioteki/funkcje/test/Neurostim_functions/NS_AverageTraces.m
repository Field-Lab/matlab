function signal=NS_AverageTraces(FileName,Timings,Channels,TimeRange,NS_GlobalConstants);
%Plots many traces showing raw data for given channel. For ewach trace the
%starting timing point is defined in 'Timings' input array.

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

%t=[TimeRange(1):TimeRange(2)]/Fs;

full_path=[pwd '\' 'data' FileName];
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
signal=zeros(TimeRange(2)-TimeRange(1)+1,length(Channels));

for i=1:length(Timings)
    data=rawFile.getData(Timings(i)+TimeRange(1),TimeRange(2)-TimeRange(1)+1);
    size(signal);
    size(data);
    %type(data)
    s=double(data(:,Channels+1));
    s1=signal+s; %first channel is a TTL channel    
    signal=s1;
end

signal=signal/length(Timings);