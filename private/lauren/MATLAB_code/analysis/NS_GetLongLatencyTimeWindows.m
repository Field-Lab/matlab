function [filtDataCat movieNumberTimes] = NS_GetLongLatencyTimeWindows(FileName, NS_GlobalConstants, windowLength, delayLength, plotLength)

% returns data as a 65 channels x nSamples matrix (first channel for
% triggers, set to zeroes)
%
%
% note that if a frequency is used in more than one movie chunk, the data from different movie chunks
% grouped into separate p-files according to the movie number, but appear in the same p-file folder

        
c = kaiserord([100 800],[0 1], [0.001 0.001], 20000, 'cell');
b = fir1(c{:});
shiftSize = floor(length(b)/2);

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
    
    timeWindows{i} = [];
    
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
    
    %[patterns,patternsIndexes,status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<NASGU>
    
    
    %% determines frequency of stimulation for this chunk
    stimTimes = MovieData(1:3:end); %every third value starting at 2
    stimPatterns = MovieData(2:3:end);
    stimInt = diff(stimTimes);
    if ~all(stimInt == stimInt(1)) || ~all(stimPatterns == stimPatterns(1))
        error('movie chunk doesn''t have constant frequency or more than one pattern is applied')
    end
    
    freq = 20000/stimInt(1);
       
    if freq == 5
        Traces=int16(zeros(nEvents, nRepeats, length(Channels), windowLength+100));
        
        for j=1:nRepeats %iterates through number of times movie is played
            
            %checks to make sure raw data file has all of the samples that the
            %stimulus files thinks it does...
            if totalSamples <= MovieBegin+repeatPeriod*j+windowLength
                warndlg(['Raw data file ' full_path ' does not have as many samples as expected by stimulus files']);
                return
            end
            
            
            % added windowLength samples to retrieved portion of data to prevent index
            % out-of-bounds error
            RawData=int16(rawFile.getData(MovieBegin + repeatPeriod*(j-1), repeatPeriod + windowLength+100)'); %extra 100 is for filtering edge effects
            
            for k = 1:nEvents
                index=(k-1)*3;
                t = MovieData(index+1);
                PatternNumber = MovieData(index+2);
                
                windowStart = MovieBegin + repeatPeriod*(j-1)+t+delayLength+1;
                windowEnd = MovieBegin + repeatPeriod*(j-1)+t + windowLength;
                timeWindows{i} = [timeWindows{i}; windowStart windowEnd]; %#ok<AGROW>
                
                Traces(k,j,:,:)=RawData(2:end, t+1:t+windowLength+100);
            end
        end
        Traces = reshape(Traces,nEvents*nRepeats,length(Channels),[]);
        meanDaq = int16(round(mean(mean(Traces(:,PatternNumber,:)))));
                
        TracesFilt = filter(b,1,double(Traces),[],3);
        
        TracesFiltInt{i} = int16(zeros(nEvents*nRepeats, length(Channels), windowLength));
        TracesFiltInt{i}(:,:,delayLength+1:windowLength) = round(TracesFilt(:,:,shiftSize+1+delayLength:shiftSize+windowLength)); %shiftSize corrects for shift generated by filtering
        
        
        if 1
            
            meanDaq = mean(mean(Traces(:,PatternNumber,:)));
            meanTrace = mean(squeeze(Traces(:,PatternNumber,:)),1);
            meanTraceFilt = mean(squeeze(TracesFiltInt{i}(:,PatternNumber,:)),1);
            
            
            figure
            hold on
            
            
            for kk = 1:nEvents*nRepeats

                
                %plot(Traces(kk,PatternNumber:)-meanDaq, 'k')
                %plot(-shiftSize:size(TracesFilt,2)-shiftSize-1, TracesFilt(kk,:), 'r')
                %plot([delayLength, delayLength], [-100 300], 'b-')
                
            end
            %plot([delayLength, delayLength], [-100 300], 'b-')
            plot([0 plotLength], [100 100], 'b-')
            plot([0 plotLength], [150 150], 'b-')
            
            %plot means for all repetitions (filtered vs unfiltered)
            plot(meanTrace(1:plotLength) - meanDaq + 100, 'k-')
            plot(meanTraceFilt(1:plotLength) + 150, 'r-')
            
            
            %plot an example trial from each
            plot(squeeze(Traces(3,PatternNumber,1:plotLength))-meanDaq, 'k-')
            plot(squeeze(TracesFiltInt{i}(3,PatternNumber,1:plotLength)), 'r-')
            
            set(gca,'ylim', [-100 300])
            
        end
        clear Traces TracesFilt
        TracesFiltInt{i} = TracesFiltInt{i}+meanDaq;
            
    end

end  

currentTime = 1;
filtDataCat = [];
movieNumberTimes = [];
for ii = 1:length(TracesFiltInt) %goes through movie numbers
    if ~isempty(TracesFiltInt{ii})
        filtDataCat = [filtDataCat reshape(shiftdim(TracesFiltInt{ii},1), 64, [])]; %#ok<AGROW>
        nOccurances = size(TracesFiltInt{ii},1);
        movieNumberTimesTemp = [ones(1,nOccurances)*ii; currentTime+(0:windowLength:windowLength*(nOccurances-1))];
        movieNumberTimes = [movieNumberTimes movieNumberTimesTemp];
        currentTime = movieNumberTimes(2,end)+windowLength;
    end
end

%trigger channel
filtDataCat = [int16(zeros(1, size(filtDataCat,2))); filtDataCat];