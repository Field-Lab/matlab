function [PDChunkNumber, MovieBegin, nRepeats, repeatPeriod]=NS_SaveRespForMovieTriplets(FileName, WritePath,...
    NS_GlobalConstants, traceLength)

% preprocesses raw electrical stimulation (STIM64 or STIM512) data
% 
% Determines the time points at which patterns are applied, then extracts a
% window of time beginning with the pattern and outputs it to a new file.
%  
% Input Arguments
% FileName - string in the form '###'
% WritePath - string, directory to write to
% NS_GlobalConstants - structure, contains SamplingFrequency, ChipAddresses, NumberOfChannelsPerChip, CurrentRanges
% traceLength - number of samples to process after the pattern is applied
%
% Output Arguments
% PDChunkNumber
% MovieBegin
% nRepeats
% repeatPeriod
%
% current directory should contain the data (down to date-piece folder)
%
% this function was originally written by Pawel, but has been modified (in mostly cosmetic ways) by
% Lauren in 2009-09
% Some comments added by Geoff Weiner 2013-01


%ArrayID=1;

% The location of the data (set to the current directory).
full_path=[pwd filesep 'data' FileName] %#ok<NOPRT>

% Creates a java object (rawFile) that let's you intereact with the raw data files (sort of the stitched together bin files).
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(full_path);

% Returns the number of samples contained in the raw data file.
totalSamples = getRawDataFileSize(full_path);

% Generates the full filenames for movie and pattern.
filename_movie=['movie' FileName] %#ok<NOPRT>
SPfilename=['pattern' FileName] %#ok<NOPRT>

% Switch to the correct path to open movie LG 1/23/2014
% tempDir = pwd; 
% cd /Volumes/Analysis/2012-09-24-3/stim_files
% keyboard; 
% Open the movie file.
fid0=fopen(filename_movie,'r','b');
% Read through the header.
header=readMHchunk(fid0);
% Number of movie chunks x number of amplitudes.
nMovies=header.number_of_movies;
%cd(tempDir); clear tempDir; 

% Calculates the number of channels, 1:64.
Channels = 1:NS_GlobalConstants.NumberOfChannelsPerChip*length(NS_GlobalConstants.ChipAddresses);
% Initialize a stopwatch.
t0 = clock;


% Iterates through all the movies.
% It always saves an extra movie at the end that should be discarded.

% Unless specific movie numbers defined!!
% Minus 1 because the last movie is expected to be empty.
for i=1:nMovies-1
    
    % Displays progress/time left.
    if i>1
        finished = (i-1)/(nMovies-1); % proportion of files created so far
        disp(sprintf('finished writing %0.1f%% of files', finished*100))
        tnow = clock;
        timeElapsed = etime(tnow, t0); % time elapsed since loop started
        estimatedTimeLeft = (timeElapsed/finished)*(1-finished);
        disp(sprintf('estimated time left: %0.1f seconds',estimatedTimeLeft))
    end
    
    % Read the chunk type ID.
    ID=fread(fid0,8,'int8')';
    
    % Checks for MD chunks by returning an error if ID is an SC chunk.
    if ID==[75 116 5 96 -84 122 -59 -64]; % SC chunk ID.
        error('command chunk found in the movie file');
        %ChunkSize = fread(fid0,1,'int64'); %read in the chunk size
        %commands=fread(fid0,ChunkSize,'int32');        
    elseif ID==[114 -69 27 -4 99 66 -12 -123] % MD chunk ID.
        
        % Advance the file pointer to the correct position for the next iteration.
        ChunkSize=fread(fid0, 1, 'int64'); % Size of MD chunk.
        fread(fid0, ChunkSize, 'int32'); % Move the pointer along that size.
        
        % Reading in the movie parameters.
        ChunkData=NS_MovieData(FileName,i,NS_GlobalConstants); %stimulus information: times at which patterns are played
        
        % Interprets first 6 values in ChunkData.
        % See NS_DecodeMovieDataChunk for descriptions.
        [PDChunkNumber,MovieBegin,nRepeats,repeatPeriod,MovieData]=NS_DecodeMovieDataChunk(ChunkData);
%         keyboard; 
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
    
    % Calculate the number of pattern applications in the current movie.
    % MovieData is structured as: [ x    y    z ...]'
    %               x = event time (in samples)
    %                   y = ID of the pattern applied
    %                       z = the number 1 (not used)
    % So every 3 entries define the information for 1 event.
    nEvents=length(MovieData)/3; 
    Events=zeros(nEvents,1);
    
    % Read the pattern information in.
    [Patterns,PatternsIndexes,status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<NASGU>
    
    % Create a status file for each movie and save it.
    FullName_status=[WritePath filesep 'status_files' filesep 'status_m' num2str(i)];
    save(FullName_status,'status')
    
    % Store the number of different patterns.MovieBegin + repeatPeriod*(j-1)
    NumberOfPatterns=length(PatternsIndexes);
    
    % Create Traces, a MxNxOxP dimensional matrix of zeros.
    %   M = nEvents
    %       N = nRepeats
    %           O = length(Channels), typcially 64
    %               P = traceLength
    Traces=int16(zeros(nEvents, nRepeats, length(Channels), traceLength));
    
    % Iterates through the repetitions of the current (i-th) movie.
    for j=1:nRepeats 
        
        % Checks to make sure raw data file has all of the samples that the stimulus files thinks it does.
        if totalSamples <= MovieBegin+repeatPeriod*j+traceLength
            warndlg(['Raw data file ' full_path ' does not have as many samples as expected by stimulus files']);
            return
        end
       
        % Gets the raw data at the time MovieBegin + j-1 repetition periods up through the length of 1 repetition plus a trace length.
        % LH: Added traceLength samples to retrieved portion of data to prevent index out-of-bounds error.
        RawData=int16(rawFile.getData(MovieBegin + repeatPeriod*(j-1), repeatPeriod + traceLength)');
        
        % Iterates through each event of the j-th repetition of the i-th movie.
        for k = 1:nEvents
            % Index moves through MovieData in steps of 3.
            index=(k-1)*3;
            t = MovieData(index+1);
            PatternNumber = MovieData(index+2);
            
            % Store the PatternNumber in a vector.
            Events(k)=PatternNumber;
            
            % Fill the Traces vector with the raw data from ALL channels.
            % First channel is the trigger channel, which is unused
            Traces(k,j,:,:)=RawData(Channels+1, t+1:t+traceLength);
        end
    end
  
    % Iterate over the patterns in the i-th movie.
    for l=1:NumberOfPatterns
        
        % Return the indices of the events in which pattern l was applied.
        WhichEvents=find(Events==l);
        
        if ~isempty(WhichEvents)
            % Extract the information about the pattern.
            Pattern=NS_ExtractPatternFromPatternDataChunk(Patterns,PatternsIndexes,l);%#ok<NASGU>
            
            % Reshape the Traces matrix by picking out the events with this pattern and collapsing the first two dimensions.
            % TracesToSave is a MxNxO matrix:
            %   M = nRepeats*length(WhichEvents)
            %       N = length(Channels), typically 64
            %           O = traceLength
            TracesToSave=reshape(Traces(WhichEvents,:,1:length(Channels),:),nRepeats*length(WhichEvents),length(Channels),traceLength);
            
            % Reshape the TracesToSave matrix so it can be written to a binary file.
            % a is a vector that contains all of the TracesToSave matrix.
            % b is a vector that contains information about the TracesToSave matrix.
            %   b(1) = nRepeats*length(WhichEvents)
            %       b(2) = length(Channels), typically 64
            %           b(3) = traceLength
            %               b(4:67) = Channels
            STTS=size(TracesToSave);
            a=reshape(TracesToSave,STTS(1)*STTS(2)*STTS(3),1);
            b=zeros(1000,1);
            b(1)=STTS(1);
            b(2)=STTS(2);
            b(3)=STTS(3);
            b(3+1:3+length(Channels))=Channels';
            
            % o is the output data that will be saved.
            o=[b' a'];

            % Make the directory for the pattern.
            if ~exist([WritePath filesep 'p' num2str(l)], 'file')
                mkdir([WritePath filesep 'p' num2str(l)])
            end

            % Write o to the p_m file.
            FullName=[WritePath filesep 'p' num2str(l) filesep 'p' num2str(l) '_m' num2str(i)];
            fid=fopen(FullName,'wb','ieee-le');
            fwrite(fid,o,'int16');
            fclose(fid);

            % Write information about the pattern to a pattern_m file.
            FullName_pattern=[WritePath filesep 'pattern_files' filesep 'pattern' num2str(l) '_m' num2str(i)];
            save(FullName_pattern,'Pattern');
            
        end
    end 
    clear Traces;
    clear TracesToSave;
end  