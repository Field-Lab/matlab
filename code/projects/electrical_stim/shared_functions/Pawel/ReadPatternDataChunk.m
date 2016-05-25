function [patterns_out,PatternsIndexes,Status]=ReadPatternDataChunk(filename,number_of_PD_chunk,NS_GlobalConstants)
%Function returns the full information about the real time data for each
%pattern in given PD chunk, and about the status of each channel of each
%chip for the time when these patterns are in use. To calculate the actual
%current values, one has to know both the real time data and the status for
%each channel.
%The output dat are:
%1) patterns_out - this is an one-dimensional array of structures. The length
%of the array is identical to total number of channels in all the patterns
%in given PD chunk. Each element of the array is a structure comprising two
%fields:
%1a) channel - channel number;
%1b) data - the array of the size 5xN. N is equal to the number of phases
%that define the pulse in this channel. The five rows include information
%about:
%- DAC values;
%- value of "record" signal;
%- value of "connect" signal;
%- value of "discharge" signal;
%- value of "hold" signal;
%2) PatternsIndexes - This is a one-dimensional array, which length is
%identical to the numbers of patterns in the PD chunk. Each value
%corresponds to one pattern and defines which stimulation pulses, defined
%in the patterns_out structure, belong to given pattern. To reconstruct the
%full n-th pattern, one must:
%2a) read the value number n from PatternsIndexes;
%2b) read in the value number n-1 from PatternsIndexes;
%2c) for each value i that is greater than value obtained in 2b, and
%less or equal to value obtained in 2a), the value number i from
%patterns_out array belongs to n-th pattern.
%Example: if the PatternIndexes array id of the form: [1 2 4 5], it means
%that:
%- there are four patterns;
%- the patterns number 1,2,4 comprise only one channel;
%- the pattern number 3 comprise 2 channels.
%
%
% Lauren's Notes
%   - only used field in Status structure is Status.ChannelStatus.range
%   (rest of values are just defaults)
%



ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

NumberOfChips=length(ChipAddresses);
NumberOfChannels=NumberOfChips*NumberOfChannelsPerChip;

ChipStatus=struct('StimActive',0','ConnectDelay',100,'RefreshDelay',100,'RecordDelay',110,'DischargeDelay',105,'HoldDelay',115);
ChipsStatus(1:NumberOfChips)=ChipStatus;

ChannelStatus=struct('active',0,'mode',0,'range',2);
ChannelsStatus(1:NumberOfChannels)=ChannelStatus;

Status=struct('ChipsStatus',ChipsStatus,'ChannelsStatus',ChannelsStatus);

fid=fopen(filename,'r','b');
header=readPHchunk(fid)


% if number_of_PD_chunk>header.number_of_chunks
%     error('the addressed number of PD chunk is too large');
% end

patterns_out=[];

PatternsIndexes=[];
pattern_index=0;

%% 

%this part seems to only serve the purpose of getting to the right place in
%the pattern file (and update Status.ChipsStatus.range values)
number=1;
disp('reading pattern data file...')
while number<number_of_PD_chunk
    ID=fread(fid,8,'int8')';
    if ID==[75 116 5 96 -84 122 -59 -64] %if this is a SC chunk...
        %mess='command chunk';
        size=fread(fid,1,'int64'); %read in the chunk size
        commands=fread(fid,size,'int32');
        for j=1:length(commands)            
            %only changes Status.ChannelStatus.range values
            Status=CommandDecoder(ChipAddresses,NumberOfChannelsPerChip,Status,commands(j));
        end
    elseif ID==[106 -23 100 -113 -38 79 0 -93]
        %mess='pattern chunk';
        size=fread(fid,1,'int64');
        %patterns=fread(fid,size,'int32');
        fread(fid,size,'int32');
        number=number+1;
    end
end

%%

ID=[1 1 1 1 1 1 1 1];
while ID~=[106 -23 100 -113 -38 79 0 -93]
    ID=fread(fid,8,'int8')';
    if ID==[75 116 5 96 -84 122 -59 -64] %if this is a SC chunk...
        size=fread(fid,1,'int64'); %read in the chunk size
        commands=fread(fid,size,'int32');
        for j=1:length(commands)                        
            Status=CommandDecoder(ChipAddresses,NumberOfChannelsPerChip,Status,commands(j));
        end
        %Status.ChannelsStatus(59)
    elseif ID==[106 -23 100 -113 -38 79 0 -93]
        size=fread(fid,1,'int64');
        patterns=fread(fid,size,'int32');
        j=1;
        
        %iterates through patterns
        patternNumber = 0;
        while j<size % for every pattern... size to ilosc wartosci (I32) w chunku a nie kanalow!!
            patternNumber = patternNumber + 1;
%             if patternNumber == 2
%                 keyboard
%             end
  
            pattern=[];
            number_of_channels=patterns(j);
            j=j+1;
            for k=1:number_of_channels % for every channel in the pattern...
%                 if (k == 27 || k == 34) && patternNumber==2
%                     keyboard
%                 end
                
                channel=patterns(j);
                j=j+1;
                number_of_phases=patterns(j);
                j=j+1;
                scaling_factor=patterns(j);
                j=j+1;
                %index=1;
                
                %durations=zeros(1,number_of_phases);
                %DAC_and_switches=zeros(1,number_of_phases);
                
                index0=1;
                ChannelData=zeros(5,number_of_phases);
                for l=1:number_of_phases % for each phase of the pulse in this channel...
                    DAC_and_switches=patterns(j);
                    j=j+1;
                    duration=patterns(j);
                    j=j+1;
                    index1=index0+duration;
                    
                    %pattern_channel(index0:index1-1)=DAC_and_switches_decode(DAC_and_switches);
                    channel_phase=DAC_and_switches_decode(DAC_and_switches);
                    ChannelData(1,index0:index1-1)=channel_phase.DAC;
                    ChannelData(2,index0:index1-1)=channel_phase.record;
                    ChannelData(3,index0:index1-1)=channel_phase.connect;
                    ChannelData(4,index0:index1-1)=channel_phase.discharge;
                    ChannelData(5,index0:index1-1)=channel_phase.hold;
                    channel_data=struct('channel',channel,'data',ChannelData);
                    
                    index0=index1;
                end
                
                pattern=[pattern channel_data]; %#ok<AGROW>
            end
            pattern_index=pattern_index+number_of_channels;
            PatternsIndexes=[PatternsIndexes pattern_index]; %#ok<AGROW>
            patterns_out=[patterns_out pattern]; %#ok<AGROW>
        end
        
    else
        error(['chunk ID:' num2str(ID) ' not correct']);
    end
end
disp('done.')
fclose(fid);