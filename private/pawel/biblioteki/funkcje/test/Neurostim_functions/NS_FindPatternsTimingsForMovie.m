function [Timings,PDChunkNumber]=NS_FindPatternsTimingsForMovie(FileName,PatternNumber,MovieNumber,NS_GlobalConstants);
%This function finds all the time points for which pulse on given electrode
%was generated (started). The function is looking for this time points in
%one movie. 
%Input arguments:
%FileName - in the form of '003';
%-Channel - number of electrode, count from 1 (the function takes care of
%the problem with numbering of channels in raw data files, where actually
%channel number one is the TTL channel; so if you want electrode number
%one, just set this variable to 1);
%-MovieNumber - which movie in given Movie File should be searched;
%-NS_GlobalConstants - the stadard structure.
%Outputs:
%-Timings: just points in time, in reference to the beginning of the run,
%not to the beginning of given movie!
%-PDChunkNumber - in which Pattern Data Chunk in corresponding Pattern File
%the patterns called by given movie are defined.

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

filename_movie=['movie' FileName];
l=length(filename_movie);
SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)];

fid0=fopen(filename_movie,'r','b');
%ID=fread(fid0,8,'int8')'
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies

offset=6; %the value 6 come sfrom the fact that 6 values on the beginning on the data array are not for events definition

%read in and ignore all the "previous" movies -this is just to find the
%position in the file where the data for THE movie begins. Should be
%changed to "skip" later.
for i=1:MovieNumber-1
    index=0;    
    ID=fread(fid0,8,'int8')';
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
        ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
        commands=fread(fid0,ChunkSize,'int32');        
    elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
        ChunkSize=fread(fid0,1,'int64');
        data=fread(fid0,ChunkSize,'int32');
    end
end

%Now, read THE movie in
ID=fread(fid0,8,'int8')';
if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
    ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
    commands=fread(fid0,ChunkSize,'int32');        
elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
    ChunkSize=fread(fid0,1,'int64');
    data=fread(fid0,ChunkSize,'int32');
    %reading in the movie parameters:
    PDChunkNumber=data(1);
    MovieBegin=data(2)*(2^30)+data(3); %beginning of first iteration of the movie
    RepetNumber=data(4); %number of repetitions of the movie;
    RepetPeriod=data(5)*(2^30)+data(6); %period of movie repetitions;   
end

PatternsNumbersOfChannel=NS_FindPatternsForChannelInChunk(FileName,PDChunkNumber,Channel,NS_GlobalConstants);
if length(PatternsNumbersOfChannel)==0
    error('No patterns for given channel found')
elseif length(PatternsNumbersOfChannel)>1
    error('More than one pattern for given channel found') 
    %however, it should be actually checked, if there is more than one pattern with this electrode sending current (and not just belonging to pattern)                                                       
end

%PulseLength=length(patterns(index).data); % length of the pulse for the given channel in sampling periods
PatternNumberOfChannel=PatternsNumbersOfChannel(1);
%PatternNumberOfChannel=min(find(PatternsIndexes>=index)); %the number of pattern where data for our channel are defined
NumberOfEvents=(ChunkSize-offset)/3; %number of events in one repetition of the movie! - all the events, not only the ones that are interesting here

TimeTriggers=[]; %the array will be loaded by the values of time, when each pulse in given
%channel starts. This 'triggers' will be then used to record adequate
%fragments of raw data to show them on the screen.

%Find all the triggers within one iteration of the movie:
for i=1:NumberOfEvents
    PatternNumber=data((i-1)*3+2+offset); 
    if PatternNumber==PatternNumberOfChannel
        Time=data((i-1)*3+1+offset)+1; %indexing in Matlab from 1!
        TimeTriggers=[TimeTriggers Time];
        ScalingFactor=data((i-1)*3+3+offset); 
        if ScalingFactor ~= 1
            warning('scaling factor defined in the movie file is different then 1!')
            disp(num2str(ScalingFactor));
        end        
    end
end 

Timings=[];
%All the triggers for all the iterations of the movie:
for i=1:RepetNumber
    Timings=[Timings TimeTriggers+MovieBegin+RepetPeriod*(i-1)];
end

fclose(fid0);