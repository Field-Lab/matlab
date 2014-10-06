function [Pulse,Status,patterns,PatternsIndexes]=NS_FindPulseShapeForMovie(FileName,Channel,PDChunkIndex,NS_GlobalConstants);
%Pulse - array of the size 5*N, where N is the length of the pulse in
%sampling periods, and the 5 rows are:
%row 1 - the DAC values;
%row 2 - values of "record" signal;
%row 3 - values of "connect" signal;
%row 4 - values of "discharge" signal;
%row 5 - values of "hold" signal.
%Status - describes the status of the chips and individual channels.
ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

SPfilename=['pattern' FileName];

[patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,PDChunkIndex,NS_GlobalConstants);
lp=length(patterns);
for i=1:lp %for each channel for which there is stimulation signal defined in this PatternDataChunk...
    if patterns(i).channel==Channel %find the index for this channel; this is index for the 'patterns' array, 
                                    %it does not define to which pattern
                                    %given channel belongs to. To find
                                    %that, one has to decode the
                                    %'PatternsIndexes'. It is done in this function several lines below. See also:
                                    %ReadPatternDataChunk help.
        if find(patterns(i).data(1,:)~=0) %if the channel is not only in this pattern, but actually sends come current...                        
            index=i;
        end                            
        %index=i;
    end
end
%PulseLength=length(patterns(index).data); % length of the pulse for the given channel in sampling periods
%PatternNumberOfChannel=min(find(PatternsIndexes>=index)); %the number of pattern where data for our channel are defined
%NumberOfEvents=(ChunkSize-offset)/3; %number of events in one repetition of the movie! - all the events, not only the ones that are interesting here

Pulse(1,:)=patterns(index).data(1,:); %this values can be potentially scaled further by defining the scaling factor parameteer value in the movie file.
Pulse(2,:)=patterns(index).data(2,:);
Pulse(3,:)=patterns(index).data(3,:);
Pulse(4,:)=patterns(index).data(4,:);
Pulse(5,:)=patterns(index).data(5,:);