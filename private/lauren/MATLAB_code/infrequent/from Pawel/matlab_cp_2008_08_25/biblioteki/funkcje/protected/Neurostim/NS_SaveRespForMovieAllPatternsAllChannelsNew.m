function [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew(FileName,WritePath,FigureProperties,NS_GlobalConstants);
% This function generates quickly plots of overlaid responses for all
% patterns and all movies FOR 1-EL. STIMULATION!!!! It is assumed that the
% number of pattern called in the movie is the number of stimulating
% electrode.
ArrayID=1;

ChipAddresses=NS_GlobalConstants.ChipAddresses;
NumberOfChannelsPerChip=NS_GlobalConstants.NumberOfChannelsPerChip;
CurrentRanges=NS_GlobalConstants.CurrentRanges;
Fs=NS_GlobalConstants.SamplingFrequency;

%full_path=[pwd '\' 'data' FileName];
full_path=[pwd filesep 'data' FileName];
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);

TraceLength=70;

filename_movie=['movie' FileName];
l=length(filename_movie);
SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)];

fid0=fopen(filename_movie,'r','b');
%ID=fread(fid0,8,'int8')'
header=readMHchunk(fid0);
NumberOfMovies=header.number_of_movies;

Channels=[1:64];
t0 = clock;

for i=1:NumberOfMovies-1 %unless specific movie numbers defined!! minus 1 because the last movie is expected to be empty
    if i>1
        finished = (i-1)/(NumberOfMovies-1); % proportion of files created so far
        disp(sprintf('finished writing %0.1f%% of files', finished*100))
        tnow = clock;
        timeElapsed = etime(tnow, t0); %time elapsed since loop started
        estimatedTimeLeft = (timeElapsed/finished)*(1-finished);
        disp(sprintf('estimated time left: %0.1f seconds',estimatedTimeLeft))
    end
        
    ID=fread(fid0,8,'int8')';
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
        error('command chunk found in the movie file');
        ChunkSize=fread(fid0,1,'int64'); %read in the chunk size
        commands=fread(fid0,ChunkSize,'int32');        
    elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
        ChunkSize=fread(fid0,1,'int64');
        ChunkData=fread(fid0,ChunkSize,'int32');
        %reading in the movie parameters:
        ChunkData=NS_MovieData(FileName,i,NS_GlobalConstants);
        [PDChunkNumber,MovieBegin,RepetNumber,RepetPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
    end
    NumberOfEvents=length(MovieData)/3
    Events=zeros(NumberOfEvents,1);
    
    [Patterns,PatternsIndexes,Status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants);
    FullName_status=[WritePath filesep 'status_m' num2str(i)];
    save(FullName_status,'Status');
            
    NumberOfPatterns=length(PatternsIndexes);
    Traces=int16(zeros(NumberOfEvents,RepetNumber,length(Channels),TraceLength));  % always 7 electrodes will be shown - should be more flexible !!!!!!
                 
    for j=1:RepetNumber
        RepIndex=j
        RawData=int16(rawFile.getData(MovieBegin+RepetPeriod*(j-1),RepetPeriod)');
        %size(RawData)
        for k=1:NumberOfEvents
            index=(k-1)*3;
            t=MovieData(index+1);
            PatternNumber=MovieData(index+2);
            Events(k)=PatternNumber;                     
            Traces(k,j,:,:)=RawData(Channels+1,t+1:t+TraceLength); %always 7 channels
        end
    end
  
    %figure(i)
    %plot(Events,'bd-');
    for l=1:NumberOfPatterns % this loop should be running only over patterns used in this movie !!!!
        WhichEvents=find(Events==l) %which events in this movie corresponded to given pattern
        if length(WhichEvents)>0
            Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,l);                                    
            TracesToSave=reshape(Traces(WhichEvents,:,1:length(Channels),:),RepetNumber*length(WhichEvents),length(Channels),TraceLength);% might be <7                 
            STTS=size(TracesToSave);        
        
            a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
            b=zeros(1000,1);
            b(1)=STTS(1);
            b(2)=STTS(2);
            b(3)=STTS(3);
            b(3+1:3+length(Channels))=Channels';
            o=[b' a'];            
            FullName=[WritePath filesep 'p' num2str(l) '_m' num2str(i)]
            fid=fopen(FullName,'wb','ieee-le');
            fwrite(fid,o,'int16');
            fclose(fid); 
            
            FullName_pattern=[WritePath filesep 'pattern' num2str(l) '_m' num2str(i)]            
            save(FullName_pattern,'Pattern');
        end
    end 
    clear Traces;
    clear TracesToSave;
end  