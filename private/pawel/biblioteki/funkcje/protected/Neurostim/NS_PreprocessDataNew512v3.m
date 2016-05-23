function [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,ClusterFileName]=NS_PreprocessDataNew512v3(FileName,WritePath,WritePathFigs,ArrayID,FigureProperties,NS_GlobalConstants,makePlots,MoviesToRead,PatternsToRead,ChannelsToRead,MovieNumberOffset);
% This function generates quickly plots of overlaid responses for all
% patterns and all movies FOR 1-EL. STIMULATION!!!! It is assumed that the
% number of pattern called in the movie is the number of stimulating
% electrode.
% PatternsToRead - if movie includes 64 patterns, but we want to process
% only some of them, the indexes should be specifiec here. The aray of
% indexes is not array of patterns really, but indexes for the array ofPatternsUsed
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

TraceLength=140;

filename_movie=['movie' FileName]
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
        i1
        PatternIndexes=find(PatternsUsed==i1); %number of events corresponding to given pattern
        PatternsIndexesToRead=[PatternsIndexesToRead' PatternIndexes']';
    end
              
    NumberOfEvents=length(PatternsIndexesToRead);  % this will not work when some patterns are used more than once in the movie!!!
    Events=zeros(NumberOfEvents,1);
    
    %c) read in corresponding pattern_data_chunk, save status into file, predefine array for RAW data  
    SPfilename
    PDChunkNumber
    NS_GlobalConstants
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants);
    'ReadPatternDataChunk done'
    FullName_status=[WritePath filesep 'status_m' num2str(i+MovieNumberOffset)];
    save(FullName_status,'Status');
            
    PatternsToRead
    NumberOfPatterns=length(PatternsIndexes)
    NumberOfPatternsMax=max(NumberOfPatternsMax,NumberOfPatterns)
    NumberOfEvents
    RepetNumber
    length(ChannelsToRead)
    TraceLength
    numberofbytes=NumberOfEvents*RepetNumber*length(ChannelsToRead)*TraceLength*2    
    if numberofbytes>0
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
    
    for l=1:length(DifferentPatterns) %1:NumberOfPatterns % this loop should be running only over patterns used in this movie !!!!
        ThisPattern=DifferentPatterns(l)        
        WhichEvents=find(Events==ThisPattern) %which events in this movie corresponded to given pattern
        NumberOfRepetitionsMax=max(NumberOfRepetitionsMax,RepetNumber*length(WhichEvents))
        if length(WhichEvents)>0            
            Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,ThisPattern);                                    
                                    
            TracesToSave=reshape(Traces(WhichEvents,:,1:length(Channels),:),RepetNumber*length(WhichEvents),length(Channels),TraceLength);              
            STTS=size(TracesToSave);        
        
            a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
            b=zeros(1000,1);
            b(1)=STTS(1);
            b(2)=STTS(2);
            b(3)=STTS(3);
            b(3+1:3+length(Channels))=Channels';
            b(4+length(Channels))=length(WhichEvents);
            o=[b' a'];            
                        
            if ArrayID==1502
                ThisPattern
                FullName=[WritePath filesep 'p' num2str(NS512_519ApplyMess(ThisPattern)) '_m' num2str(i+MovieNumberOffset)];
            else
                FullName=[WritePath filesep 'p' num2str(ThisPattern) '_m' num2str(i+MovieNumberOffset)];
            end
            fid=fopen(FullName,'wb','ieee-le');
                                    
            fwrite(fid,o,'int16');
            fclose(fid);                    
                        
            if ArrayID==1502
                FullName_pattern=[WritePath filesep 'pattern' num2str(NS512_519ApplyMess(ThisPattern)) '_m' num2str(i+MovieNumberOffset)];
            else
                FullName_pattern=[WritePath filesep 'pattern' num2str(ThisPattern) '_m' num2str(i+MovieNumberOffset)];  
            end                                   
            save(FullName_pattern,'Pattern');
            
            %part below - only for 1-el. scan!!
            if makePlots==1
                ChannelsToPlot=electrodeMap.getAdjacentsTo(ThisPattern,1);                
                TracesToShow=TracesToSave(:,ChannelsToPlot,:);
                %c=int16(mean(TracesToShow));
                %for ci=1:STTS(1)
                %    TracesToShow(ci,:,:)=TracesToShow(ci,:,:)-c;
                %end
            
                y=NS_PlotManySignaturesOnArrayLayoutNew(TracesToShow,ChannelsToPlot,ones(1,STTS(1)),ArrayID,FigureProperties,NS_GlobalConstants);
                hj=gcf;
                %set(hj, 'PaperOrientation', 'portrait');
                              
                set(hj,'PaperUnits','inches');
                set(hj,'PaperSize',[13 9]);
                set(hj,'PaperPosition',[0 0 13 9]);
                %print(hc, '-dtiff', '-r120', name);        
                
                FullName=[WritePathFigs '\' 'p' num2str(ThisPattern) '_m' num2str(i)];
                print(hj, '-dtiff', '-r120', FullName);
            end
        end
    end                                                                          
    clear Traces;
    clear TracesToSave;
    end
end  
fclose(fid0);
%create cluster file:
ClusterFileName=NS_CreateClusterFile(WritePath,FileName,NumberOfMovies,NumberOfPatternsMax,NumberOfRepetitionsMax);