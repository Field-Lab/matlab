function [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,makePlots,MoviesToRead,PatternsToRead,ChannelsToRead);
% This function generates quickly plots of overlaid responses for all
% patterns and all movies FOR 1-EL. STIMULATION!!!! It is assumed that the
% number of pattern called in the movie is the number of stimulating
% electrode.
% PatternsToRead - if movie includes 64 patterns, but we want to process
% only some of them, the indexes should be specifiec here. The aray of
% indexes is not array of patterns really, but indexes for the array of
% patterns that will be loaded by this function from the stimulus files.

%ArrayID=1;

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

%full_path=[pwd '\' 'data' FileName];
full_path=[pwd filesep 'data' FileName]
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(ArrayID);

TraceLength=70;

filename_movie=['movie' FileName];
l=length(filename_movie);
SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)];

fid0=fopen(filename_movie,'r','b');
%ID=fread(fid0,8,'int8')'
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies

Channels=ChannelsToRead;
t0 = clock;
NumberOfPatternsMax=0;
NumberOfRepetitionsMax=0;

for i=MoviesToRead %1:8:NumberOfMovies-1 %unless specific movie numbers defined!! minus 1 because the last movie is expected to be empty. This loop repeats once for each amplitude (in case of scan data).    
    %a ) estimate how much time is left to complete the loop
    if i>1 
        finished = (i-1)/(NumberOfMovies-1); % proportion of files created so far
        disp(sprintf('finished writing %0.1f%% of files', finished*100))
        tnow = clock;
        timeElapsed = etime(tnow, t0); %time elapsed since loop started
        estimatedTimeLeft = (timeElapsed/finished)*(1-finished);
        disp(sprintf('estimated time left: %0.1f seconds',estimatedTimeLeft))
    end

    %b) read in single movie_data_chunk    
    ID=fread(fid0,8,'int8')';
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC (command) chunk...
        error('command chunk found in the movie file');
        ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
        commands=fread(fid0,ChunkSize,'int32');        
    elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD (movie data) chunk
        ChunkSize=fread(fid0,1,'int64');
        ChunkData=fread(fid0,ChunkSize,'int32');
        %reading in the movie parameters:
        ChunkData=NS_MovieData(FileName,i,NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);        
    end    

    %calculate number of events related to specified patterns
    PatternsUsed=MovieData(2:3:length(MovieData))
    
    PatternsIndexesToRead=[];
    for i1=PatternsToRead
        PatternIndexes=find(PatternsUsed==i1); %number of events corresponding to given pattern
        PatternsIndexesToRead=[PatternsIndexesToRead' PatternIndexes']';
    end
              
    NumberOfEvents=length(PatternsIndexesToRead)  % this will not work when some patterns are used more than once in the movie!!!
    Events=zeros(NumberOfEvents,1);
    
    %c) read in corresponding pattern_data_chunk, save status into file, predefine array for RAW data    
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants);
    'ReadPatternDataChunk done'
    FullName_status=[WritePath filesep 'status_m' num2str(i)];
    save(FullName_status,'Status');
            
    NumberOfPatterns=length(PatternsIndexes)
    NumberOfPatternsMax=max(NumberOfPatternsMax,NumberOfPatterns);
    numberofbytes=NumberOfEvents*RepetNumber*length(ChannelsToRead)*TraceLength*2
    
    Traces=int16(zeros(NumberOfEvents,RepetNumber,length(ChannelsToRead),TraceLength));                    
    
    %d) read in RAW data for each repetition, for each stimulation event
    %identify pattern number, and save RAW data coresponding to each event into
    %Traces array. We assume each repetition of the movie includes identical
    %collection of patterns (logical, isn't it).

    for j=1:RepetNumber % one iteration of this loop takes less than 0.5s for 64-channel data (122 events)                
        RepIndex=j
        RawData=int16(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod)');                        
        for k=1:NumberOfEvents % later on, single iteration should read and organize data for all ... of given pattern in this movie at once (PH, 2010-08-26)
            PatternNumber=PatternsUsed(PatternsIndexesToRead(k));
            t=MovieData(1+(PatternsIndexesToRead(k)-1)*3);                                                                
            Events(k)=PatternNumber;            
            Traces(k,j,:,:)=RawData(Channels+1,t+1:t+TraceLength); %this takes less than 0.1ms for 64-channel data                                 
        end        
    end    
    
    DifferentPatterns=unique(Events);
        
    clear RawData;
    Events
    for kl=1:NumberOfEvents                
        ThisPattern=Events(kl);
        Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,ThisPattern);   
        
        TracesToSave=reshape(Traces(kl,:,1:length(Channels),:),RepetNumber,length(Channels),TraceLength); 
        STTS=size(TracesToSave);
        a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
            b=zeros(1000,1);
            b(1)=STTS(1);
            b(2)=STTS(2);
            b(3)=STTS(3);
            b(3+1:3+length(Channels))=Channels';
            b(4+length(Channels))=1;
            o=[b' a'];   
        
        
        hj=find(Events==ThisPattern);
        PulseRep=kl-hj(1)+1;
        FullName=[WritePath filesep 'p' num2str(ThisPattern*100+PulseRep) '_m' num2str(i)];
        
        fid=fopen(FullName,'wb','ieee-le');                                    
        fwrite(fid,o,'int16');
        fclose(fid);   
        
        FullName_pattern=[WritePath filesep 'pattern' num2str(ThisPattern*100+PulseRep) '_m' num2str(i)];            
        save(FullName_pattern,'Pattern');
    end                
    
    %end                                                                          
    clear Traces;
    clear TracesToSave;
end  
fclose(fid0);
%create cluster file:
ClusterFileName=NS_CreateClusterFile(WritePath,FileName,NumberOfMovies,NumberOfPatternsMax,NumberOfRepetitionsMax);