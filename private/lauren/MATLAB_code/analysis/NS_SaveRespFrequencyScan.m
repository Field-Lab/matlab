function [PDChunkNumber, MovieBegin, nRepeats, repeatPeriod]=NS_SaveRespFrequencyScan(FileName, WritePath, NS_GlobalConstants, traceLength)

% preprocesses raw electrical stimulation (STIM64), partitioning traces the same way as the standard
% preprocessing, but giving traces applied at different frequencies different 'p-file' names
%
% ***this only works if the stimulus is structured so that each movie chunk corresponds to a
% different frequency!
%
% note that if a frequency is used in more than one movie chunk, the data from different movie chunks
% grouped into separate p-files according to the movie number, but appear in the same p-file folder
%
% function determines frequency of pulses and checks that all pulses in the movie are the same
%
% ordering of traces: goes through all pulses in sequence (in order) of first
% movie repetition, then all pulses in sequence of second movie
% repetitions, etc.
%
% this function is a modified verson of a function originally written by Pawel
% author: Lauren
%
% last modification 2010-12-28
% verified that frequencies identified correctly by setting traceLength to 300




full_path=[pwd filesep 'data' FileName] %#ok<NOPRT>
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

totalSamples = getRawDataFileSize(full_path);


filename_movie=['movie' FileName] %#ok<NOPRT>
l=length(filename_movie);
SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)] %#ok<NOPRT>

fid0=fopen(filename_movie,'r','b');
header=readMHchunk(fid0);
nMovies=header.number_of_movies;

Channels=1:64;
t0 = clock;

% iterates through movies (movies chunks x stimulus amplitudes)
for i=1:nMovies-1 %unless specific movie numbers defined!! minus 1 because the last movie is expected to be empty
    
    %displays progress/time left
    if i>1
        finished = (i-1)/(nMovies-1); % proportion of files created so far
        disp(sprintf('finished writing %0.1f%% of files', finished*100))
        tnow = clock;
        timeElapsed = etime(tnow, t0); %time elapsed since loop started
        estimatedTimeLeft = (timeElapsed/finished)*(1-finished);
        disp(sprintf('estimated time left: %0.1f seconds',estimatedTimeLeft))
    end
        
    ID=fread(fid0,8,'int8')';
    
    if ID==[75 116 5 96 -84 122 -59 -64]; %if this is a SC chunk...
        error('command chunk found in the movie file');
        %ChunkSize = fread(fid0,1,'int64'); %read in the chunk size
        %commands=fread(fid0,ChunkSize,'int32');        
    elseif ID==[114 -69 27 -4 99 66 -12 -123] %MD chunk
        
        % need to leave these in so that the file pointer is in the correct position for the next
        % iteration
        ChunkSize=fread(fid0, 1, 'int64'); %size of MD chunk
        fread(fid0, ChunkSize, 'int32'); %just to move file pointer along
        
        %reading in the movie parameters:
        ChunkData=NS_MovieData(FileName,i,NS_GlobalConstants); %stimulus information: times at which patterns are played
        
        % interprets first 6 values in chunkdata
        % see NS_DecodeMovieDataChunk for descriptions
        [PDChunkNumber,MovieBegin,nRepeats,repeatPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);

        
    end
    nEvents=length(MovieData)/3; %number of pattern applications in movie
    Events=zeros(nEvents,1);
    
    [patterns,patternsIndexes,status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<NASGU>
    FullName_status=[WritePath filesep 'status_files' filesep 'status_m' num2str(i)];
    save(FullName_status,'status')
    
    
    %% determines frequency of stimulation for this chunk
    stimTimes = MovieData(1:3:end); %every third value starting at 2
    stimPatterns = MovieData(2:3:end);
    stimInt = diff(stimTimes);
    if ~all(stimInt == stimInt(1)) || ~all(stimPatterns == stimPatterns(1))
        error('movie chunk doesn''t have constant frequency or more than one pattern is applied')
    end
    
    freq = 20000/stimInt(1);
    
    
    
    NumberOfPatterns=length(patternsIndexes);
    Traces=int16(zeros(nEvents, nRepeats, length(Channels), traceLength));
        
    for j=1:nRepeats %iterates through number of times movie is played
        
        %checks to make sure raw data file has all of the samples that the
        %stimulus files thinks it does...
        if totalSamples <= MovieBegin+repeatPeriod*j+traceLength
            warndlg(['Raw data file ' full_path ' does not have as many samples as expected by stimulus files']);
            return
        end
       
        
        % added traceLength samples to retrieved portion of data to prevent index
        % out-of-bounds error
        RawData=int16(rawFile.getData(MovieBegin + repeatPeriod*(j-1), repeatPeriod + traceLength)');
        
        for k = 1:nEvents
            index=(k-1)*3;
            t = MovieData(index+1);
            PatternNumber = MovieData(index+2);
            
            
            Events(k)=PatternNumber;
            Traces(k,j,:,:)=RawData(Channels+1, t+1:t+traceLength);
        end
    end
  
    %figure(i)
    %plot(Events,'bd-');
    for l=1:NumberOfPatterns % this loop should be running only over patterns used in this movie !!!!
        WhichEvents=find(Events==l); %which events in this movie corresponded to given pattern (should be all of them)
        if ~isempty(WhichEvents)
   
            Pattern=NS_ExtractPatternFromPatternDataChunk(patterns,patternsIndexes,l);%#ok<NASGU>

            TracesToSave=reshape(Traces(WhichEvents,:,1:length(Channels),:),nRepeats*length(WhichEvents),length(Channels),traceLength);
            STTS=size(TracesToSave);

            a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
            b=zeros(1000,1);
            b(1)=STTS(1);
            b(2)=STTS(2);
            b(3)=STTS(3);
            b(3+1:3+length(Channels))=Channels';
            o=[b' a'];

            if ~exist([WritePath filesep 'p' num2str(l) '_f' num2str(freq)], 'file')
                mkdir([WritePath filesep 'p' num2str(l) '_f' num2str(freq)])
            end

            FullName=[WritePath filesep 'p' num2str(l) '_f' num2str(freq) filesep...
                'p' num2str(l) '_f' num2str(freq) '_m' num2str(i)];
            fid=fopen(FullName,'wb','ieee-le');
            fwrite(fid,o,'int16');
            fclose(fid);

            FullName_pattern=[WritePath filesep 'pattern_files' filesep 'pattern' num2str(l) '_f' num2str(freq) '_m' num2str(i)];
            save(FullName_pattern,'Pattern');
            
        end
    end 
    clear Traces;
    clear TracesToSave;
end  