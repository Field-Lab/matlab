

NS_GlobalConstants = struct('SamplingFrequency',20000,'ChipAddresses',[30 31],'NumberOfChannelsPerChip',32,'CurrentRanges',[0.066 0.266 1.07 4.25 16.9 67.1 264 1040]);

cd /Volumes/Ripple/Data/2010-11-22-5/

FileName = '008';
    


filename_movie=['movie' FileName] %#ok<NOPRT>
l=length(filename_movie);
SPfilename=[filename_movie(1:l-8) 'pattern' filename_movie(l-2:l)] %#ok<NOPRT>

fid0=fopen(filename_movie,'r','b');
header=readMHchunk(fid0);
nMovies=header.number_of_movies;

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
        %[PDChunkNumber,MovieBegin,nRepeats,repeatPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);

        
    end
    
    %if i == 55;
        
    
    
        nEvents=length(MovieData)/3; %number of pattern applications in movie
        Events=zeros(nEvents,1);
        
        [Patterns,PatternsIndexes,status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<NASGU>
        
        
        NumberOfPatterns=length(PatternsIndexes);
        
        
%         for j=1:nRepeats %iterates through number of times movie is played
%             for k = 1:nEvents
%                 index=(k-1)*3;
%                 PatternNumber = MovieData(index+2);
%                 Events(k)=PatternNumber;
%             end
%         end
        
        for l=1:NumberOfPatterns % this loop should be running only over patterns used in this movie !!!!
            %if ~isempty(find(Events==l, 1)); %which events in this movie corresponded to given pattern                
                Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,l);%#ok<NASGU>
                
                if all(Pattern(:,1)==0) && size(Pattern,1) == 9
                    keyboard
                end
            %end
        end
    %end

end  