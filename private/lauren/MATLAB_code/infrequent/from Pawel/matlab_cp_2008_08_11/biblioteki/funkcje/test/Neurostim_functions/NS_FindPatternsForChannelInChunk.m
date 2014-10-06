function PatternsNumbersOfChannel=NS_FindPatternsForChannelInChunk(FileName,PDChunkNumber,Channel,NS_GlobalConstants);
%This functions finds numbers of all the patterns in the given pattern data
%chunk, that include given Channel. Input data:
%FileName - only the number, example: '003';
%PDChunkNumber - number of Patter Data Chunk in given file;
%Channel - number of channel.
%In principles there may be more than one pattern that include given
%channel in the same Pattern Data Chunk.

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

SPfilename=['pattern' num2str(FileName)]

[patterns,PatternsIndexes,status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants);
lp=length(patterns);
indexes=[];
for i=1:lp %for each channel for which there is stimulation signal defined in this PatternDataChunk...
    if patterns(i).channel==Channel %find the index for this channel; this is index for the 'patterns' array, 
                                    %it does not define to which pattern
                                    %given channel belongs to. To find
                                    %that, one has to decode the
                                    %'PatternsIndexes'. It is done in this function several lines below. See also:
                                    %ReadPatternDataChunk help.
        if find(patterns(i).data(1,:)~=0) %if the channel is not only in this pattern, but actually sends come current...            
            indexes=[indexes i];
        end
    end
end
indexes
%PulseLength=length(patterns(index).data); % length of the pulse for the given channel in sampling periods
for i=1:length(indexes)
    PatternsNumbersOfChannel(i)=min(find(PatternsIndexes>=indexes(i))); %the number of pattern where data for our channel are defined
end