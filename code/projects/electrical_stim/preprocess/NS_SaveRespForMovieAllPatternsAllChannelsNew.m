function [PDChunkNumber, movieBegin, numRepetitions, repetitionPeriod]= ...
    NS_SaveRespForMovieAllPatternsAllChannelsNew(rawDataPath, WritePath,...
    NS_GlobalConstants, traceLength)
%
% This function preprocesses raw electrical stimulation data
% 
% NS_SaveRespForMovieAllPatternsAllChannelsNew(...) will determine the 
% time points at which stimulation patterns are applied, extract the
% chunk of raw data corresponding to a pattern, and save the chunk to a 
% new file.
%  
% Input Arguments
% rawDataPath - full path to the filenames string in the form '###'
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
% this function was originally written by Pawel, but has been modified (in mostly cosmetic ways) by
% Lauren Jepson in 2009-09
% Some comments added by Geoff Weiner 2013-01
% Lauren Grosberg modified to read the pattern and movie files from 
% arbitrary directories 2014-10


% Creates a java object (rawFile) that let's you intereact with the raw data files (sort of the stitched together bin files).
rawFile=edu.ucsc.neurobiology.vision.io.RawDataFile(rawDataPath);
header = rawFile.getHeader();
% Returns the number of samples contained in the raw data file.
totalSamples = header.getNumberOfSamples(); 

[parentstr,datarun,~] = fileparts(rawDataPath); 
datarunNumber = datarun(5:end); 
altFileOrg = fullfile(parentstr,'Electrical', 'Output');
if ~isempty(dir(fullfile(parentstr,'movie*')))
    movieFileName= fullfile(parentstr,['movie' datarunNumber]);
    SPfilename    = fullfile(parentstr,['pattern' datarunNumber]);
elseif ~isempty(dir(fullfile(altFileOrg,'movie*')))
    movieFileName= fullfile(altFileOrg, ['movie' datarunNumber]);
    SPfilename    = fullfile(altFileOrg,['pattern' datarunNumber]);
else
    % User selects directory with the labview output files
    correctFileOrg = uigetdir(parentstr, 'Select directory containing LabVIEW output files ''movie*'' and ''pattern*'''); 
    if correctFileOrg
        movieFileName = fullfile(correctFileOrg,['movie' datarunNumber]);
        SPfilename     = fullfile(correctFileOrg,['pattern' datarunNumber]);
    else
        err = MException('MATLAB:missingLabVIEWFiles', ...
            'Cannot find labview output files that define the times and patterns of electrical stimulation');
        throw(err);
    end
end

% Open the movie file.
fid0=fopen(movieFileName,'r','b'); 
% Read through the header.
header=readMHchunk(fid0);
% Number of movie chunks x number of amplitudes.
nMovies=header.number_of_movies;

% Calculates the number of channels
Channels = 1:NS_GlobalConstants.NumberOfChannelsPerChip*length(NS_GlobalConstants.ChipAddresses);
% Initialize a stopwatch.
t0 = tic;

for i=1:nMovies-1 % the last movie is expected to be empty.
    
    % Displays progress/time left.
    if i>1
        finished = (i-1)/(nMovies-1); % proportion of files created so far
        fprintf('finished writing %0.1f%% of files', finished*100)
        estimatedTimeLeft = (toc(t0)/finished)*(1-finished);
        fprintf('estimated time left: %0.1f minutes',estimatedTimeLeft/60)
    end
    
    % Read the chunk type ID.
    ID=fread(fid0,8,'int8')';
    
    % Checks for MD chunks by returning an error if ID is an SC chunk.
    if ID==[75 116 5 96 -84 122 -59 -64]; % SC chunk ID.
        error('command chunk found in the movie file');
        %ChunkSize = fread(fid0,1,'int64'); %read in the chunk size
        %commands=fread(fid0,ChunkSize,'int32');        
    elseif ID==[114 -69 27 -4 99 66 -12 -123] % movie data chunk ID.
        % Advance the file pointer to the correct position for the next iteration.
        ChunkSize=fread(fid0, 1, 'int64'); % Size of movie data chunk.
        fread(fid0, ChunkSize, 'int32'); % Move the pointer along that size.
        % Reading in the movie parameters.
        %stimulus information: times at which patterns are played
        [PDChunkNumber,movieBegin,numRepetitions,repetitionPeriod,data] = getMovieData(movieFileName, i); 
    end
    
    % Calculate the number of pattern applications in the current movie.
    % MovieData is structured as: [ x    y    z ...]'
    %               x = event time (in samples)
    %                   y = ID of the pattern applied
    %                       z = the number 1 (not used)
    % So every 3 entries define the information for 1 event.
    nEvents=length(data)/3; 
    Events=zeros(nEvents,1);
    
    % Read the pattern information in.
    [Patterns,PatternsIndexes,status]=ReadPatternDataChunk(SPfilename,PDChunkNumber,NS_GlobalConstants); %#ok<NASGU> 
    
    % Create a status file for each movie and save it.
    FullName_status= fullfile(WritePath,'status_files',['status_m' num2str(i)]);
    save(FullName_status,'status')
    
    % Store the number of different patterns.MovieBegin + repeatPeriod*(j-1)
    NumberOfPatterns=length(PatternsIndexes);
    
    % Create Traces, a MxNxOxP dimensional matrix of zeros.
    %   M = nEvents
    %       N = nRepeats
    %           O = length(Channels), typcially 64
    %               P = traceLength
    Traces=int16(zeros(nEvents, numRepetitions, length(Channels), traceLength));
    
    % Iterates through the repetitions of the current (i-th) movie.
    for j=1:numRepetitions 
        
        % Checks to make sure raw data file has all of the samples that the stimulus files thinks it does.
        if totalSamples <= movieBegin+repetitionPeriod*j+traceLength
            err = MException('MATLAB:rawDataMissingSamples', ...
                ['Raw data file ' rawDataPath ' does not have as many samples as expected by stimulus files']);
            throw(err);
        end
       
        % Gets the raw data at the time MovieBegin + j-1 repetition periods up through the length of 1 repetition plus a trace length.
        % LH: Added traceLength samples to retrieved portion of data to prevent index out-of-bounds error.
        RawData=int16(rawFile.getData(movieBegin + repetitionPeriod*(j-1), repetitionPeriod + traceLength)');
        
        % Iterates through each event of the j-th repetition of the i-th movie.
        for k = 1:nEvents
            % Index moves through MovieData in steps of 3.
            index=(k-1)*3;
            t = data(index+1);
            PatternNumber = data(index+2);
            
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
            TracesToSave=reshape(Traces(WhichEvents,:,1:length(Channels),:),numRepetitions*length(WhichEvents),length(Channels),traceLength);
            
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