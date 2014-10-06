function Traces=NS_ReadManyTracesFromRaw(FileName,Channels,Timings,TimeStart,NumberOfSamples,Offsets,NS_GlobalConstants);
%Reads from the RAW data file number of traces from multpile channels.
%Input:
%TimeStart - what should be the number of first sample for each traces in
%reference to value of Timings for this trace
%Output:
%Traces - array of Size CxNxL, where: C - number of channels, N - number of
%traces for each channel, L - number of samples per each trace.
if min(Timings)<-TimeStart
    error('Value of TimeStart is too negative');
end

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

Traces=zeros(length(Timings),length(Channels),NumberOfSamples);
size(Traces);

full_path=[pwd '\' 'data' FileName];
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

for i=1:length(Timings)
    data=double(rawFile.getData(Timings(i)+TimeStart,NumberOfSamples)');
    %size(data)
    Traces(i,:,:)=data(Channels+1,:); %first channel is a TTL channel            
end

for i=1:length(Channels)
    Traces(:,i,:)=Traces(:,i,:)-Offsets(i);
end