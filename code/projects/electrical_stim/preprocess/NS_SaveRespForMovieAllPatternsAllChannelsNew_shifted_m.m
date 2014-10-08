function [PDChunkNumber, MovieBegin, nRepeats, repeatPeriod]=NS_SaveRespForMovieAllPatternsAllChannelsNew_shifted_m(FileName, WritePath, NS_GlobalConstants, traceLength, mShift)

% preprocesses raw electrical stimulation (STIM64) data, assigning movie
% numbers that are offset by the value mShift (useful when one dataset is
% collected in two chunks, with the second started where the previous left off)
%
% this function was originally written by Pawel, but has been modified (in mostly cosmetic ways) by
% Lauren in 2011-06
%
% 

%ArrayID=1;

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
    mForSaving = i + mShift;
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

        
        %% temporary hack!!!
%         repeatPeriod = 40000;
%         if i==1
%             MovieBegin = 440000;
%         elseif i==2
%             MovieBegin = 1440000;    
%         else
%             MovieBegin = repeatPeriod*(25*i-15);
%         end
    end
    nEvents=length(MovieData)/3; %number of pattern applications in movie
    Events=zeros(nEvents,1);
    
    [Patterns,PatternsIndexes,status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<NASGU>
    FullName_status=[WritePath filesep 'status_files' filesep 'status_m' num2str(mForSaving)];
    save(FullName_status,'status')
    
    
    NumberOfPatterns=length(PatternsIndexes);
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
        WhichEvents=find(Events==l); %which events in this movie corresponded to given pattern
        if ~isempty(WhichEvents)
             
            pForSaving = l;
 
            if pForSaving ~= 0
                
                Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,l);%#ok<NASGU>
                
                TracesToSave=reshape(Traces(WhichEvents,:,1:length(Channels),:),nRepeats*length(WhichEvents),length(Channels),traceLength);
                STTS=size(TracesToSave);

                a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
                b=zeros(1000,1);
                b(1)=STTS(1);
                b(2)=STTS(2);
                b(3)=STTS(3);
                b(3+1:3+length(Channels))=Channels';
                o=[b' a'];

                if ~exist([WritePath filesep 'p' num2str(pForSaving)], 'file')
                    mkdir([WritePath filesep 'p' num2str(pForSaving)])
                end

                FullName=[WritePath filesep 'p' num2str(pForSaving) filesep 'p' num2str(pForSaving) '_m' num2str(mForSaving)];
                fid=fopen(FullName,'wb','ieee-le');
                fwrite(fid,o,'int16');
                fclose(fid);

                FullName_pattern=[WritePath filesep 'pattern_files' filesep 'pattern' num2str(pForSaving) '_m' num2str(mForSaving)];
                save(FullName_pattern,'Pattern');
            end
            
        end
    end 
    clear Traces;
    clear TracesToSave;
end  