function [PDChunkNumber, MovieBegin, nRepeats, repeatPeriod]=NS_SaveRespEachOccurrenceSeparately(FileName, WritePath, NS_GlobalConstants, traceLength)

% preprocesses raw electrical stimulation (STIM64) data, saving traces associated with each
% occurance of a pattern within the movie file (across all movie chunks) into a separate 'p-file'
%
% function determines when new movie corresponds with another movie chunk at the same master
% amplitude or the first movie chunk at a new master amplitude, and groups across amplitudes, but
% not across movie chunks, into p-files
%
% this determination is made by comparing the current movie's at the movie data (timing of each
% pattern application) to the first movie's movie data, and checking to make sure the pulse
% amplitude of the first pattern in the movie has changed (new master amplitude)
%
% this function is a modified verson of a function originally written by Pawel
% author: Lauren
%
% last modification 2010-12-27
% needs to be more thoroughly tested before it is trusted!
%
%
%



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

occurrenceCount = 0;

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
    
    
    %determines whether this movie corresponds with the beginning the set of movie chunks at a new amplitude
    %checks if same moviedata has been observed before (same pattern application times)
        
    %get pulse amplitude of first pattern used in this movie, which is second element in MovieData (arbitrary but consistent)
    
    thisPattern = NS_ExtractPatternFromPatternDataChunk(patterns,patternsIndexes,MovieData(2));%#ok<NASGU>
        
    %find first electrode that has stimulation
    stimAmp = 0;
    for pp = 1:length(thisPattern)
        if ~all(thisPattern(pp).data(1,:) == 0) %if stimulation on this electrode
            currentStep = NS_GlobalConstants.CurrentRanges(status.ChannelsStatus(thisPattern(pp).channel).range + 1)/127;
            
            %determine pulse amplitude on this electrode in this pattern
            for jj = 1:size(thisPattern(pp).data,2)
                stimAmp = max(stimAmp, abs(thisPattern(pp).data(1,jj)*thisPattern(pp).data(3,jj)*currentStep));
            end
            break
        end
    end
    
    if stimAmp == 0
        error('problem determining whether current movie is start of new amplitude (first pattern in PD chunk has no stim)')
    end
    
    %newAmplitude = false;
    if i == 1
        firstMovieData = MovieData;
        prevStimAmp = stimAmp;
    elseif isequal(MovieData, firstMovieData)
        if prevStimAmp ~= stimAmp
            %newAmplitude = true;
            occurrenceCount = 0; %resets count of pattern occurances for new pulse amplitude
            prevStimAmp = stimAmp;
        else
            disp('warning -- first movie chunk appears to be repeated at same current value: FATAL ERROR')
            continue
        end
    end
   
    

    NumberOfPatterns=length(patternsIndexes); %note: this may include patterns that were not used in this movie
    Traces=zeros(nEvents, nRepeats, length(Channels), traceLength, 'int16');
    
    for j = 1:nRepeats %iterates through number of times movie is played
        
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
            
            
            Events(k)=PatternNumber; %keeps track of which stimulus applications correspond with which patterns
            Traces(k,j,:,:)=RawData(Channels+1, t+1:t+traceLength);
        end
    end
  
    %loop through patterns
    for l = 1:NumberOfPatterns % this loop should be running only over patterns used in this movie !!!!
        WhichEvents=find(Events==l); %which events in this movie corresponded to given pattern
        if ~isempty(WhichEvents)
            pattern = NS_ExtractPatternFromPatternDataChunk(patterns, patternsIndexes, l);%#ok<NASGU>
            
            for mm = 1:length(WhichEvents) %save each pattern occurance separately
                
                if length(occurrenceCount) < l
                    occurrenceCount(l) = 1;
                else
                    occurrenceCount(l) = occurrenceCount(l)+1;
                end
                
                TracesToSave=reshape(Traces(WhichEvents(mm),:,1:length(Channels),:),nRepeats,length(Channels),traceLength);
                
                STTS=size(TracesToSave);
                a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
                b=zeros(1000,1);
                b(1)=STTS(1);
                b(2)=STTS(2);
                b(3)=STTS(3);
                b(3+1:3+length(Channels))=Channels';
                o=[b' a'];
                
                if ~exist([WritePath filesep 'p' num2str(l) '_o' num2str(occurrenceCount(l))], 'file')
                    mkdir([WritePath filesep 'p' num2str(l) '_o' num2str(occurrenceCount(l))])
                end
                
                FullName=[WritePath filesep 'p' num2str(l) '_o' num2str(occurrenceCount(l))... 
                    filesep 'p' num2str(l) '_o' num2str(occurrenceCount(l)) '_m' num2str(i)];
                fid=fopen(FullName,'wb','ieee-le');
                fwrite(fid,o,'int16');
                fclose(fid);
                
                % write pattern file (yes they are identical for each pattern occurance)
                FullName_pattern=[WritePath filesep 'pattern_files' filesep 'pattern'...
                    num2str(l) '_o' num2str(occurrenceCount(l)) '_m' num2str(i)];
                save(FullName_pattern,'pattern');
            end
            
        end
    end 
    clear Traces;
    clear TracesToSave;
end  