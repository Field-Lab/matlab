%% Script to view movies from 2014-08-20-1/data002 

system        = 'stim512'; 
traceLength   = 100;
pathToLabViewOutputFiles = '/Volumes/Data/2014-08-20-1/'; 
datarun                  = '002';

originalDataFilePath= ['/Volumes/Data/2014-08-20-1/data' datarun];
originalRawFile     = edu.ucsc.neurobiology.vision.io.RawDataFile(originalDataFilePath);
refHeader           = originalRawFile.getHeader();
totalSamples        = refHeader.getNumberOfSamples;
Channels            = 1:(refHeader.getNumberOfElectrodes-1); %-1 because first channel has TTL pulses

if strcmpi(system, 'stim512')
    ChipAddresses=31:-1:24;
    NumberOfChannelsPerChip=64;
    ArrayID=500;
elseif strcmpi(system, 'stim64')
    ChipAddresses=[30 31];
    NumberOfChannelsPerChip=32;
    ArrayID=1;
end
% Upper limits of current ranges. Constant within a movie chunk for each electrode)
CurrentRanges=[0.066 0.266 1.07 4.25 16.9 67.1 264 1040]; % units are in uA, i.e., 66nA, 266nA, 1.07uA, 4.25uA,....
Fs = 20000; % sampling rate (in Hz)
NS_GlobalConstants=struct('SamplingFrequency',Fs,'ChipAddresses',ChipAddresses,...
    'NumberOfChannelsPerChip',NumberOfChannelsPerChip,'CurrentRanges',CurrentRanges);

% Generates the full filenames for movie and pattern output files.
fnameMovie      = [pathToLabViewOutputFiles 'movie' datarun]; 
fnamePattern    = [pathToLabViewOutputFiles 'pattern' datarun]; 

fid0    = fopen(fnameMovie,'r','b');      % Open the movie file.
header  = readMHchunk(fid0);              % Read through the header.
nMovies = header.number_of_movies;        % # movie chunks x # amplitudes

cd(pathToLabViewOutputFiles);


for i=1:nMovies-1 % Last movie is expected to be empty.
    ID  = fread(fid0,8,'int8')'; % Read the chunk type ID.
    % Checks for MD chunks by returning an error if ID is an SC chunk.
    if ID==[75 116 5 96 -84 122 -59 -64];            % SC chunk ID.
        error('command chunk found in the movie file');   
    elseif ID==[114 -69 27 -4 99 66 -12 -123]  % MD chunk ID.
        % Advance the file pointer to the correct position for the next iteration.
        ChunkSize  = fread(fid0, 1, 'int64'); % Size of MD chunk.
        fread(fid0, ChunkSize, 'int32'); % Move the pointer along that size.
        % Read in the movie parameters for both runs.
        ChunkData=NS_MovieData(datarun,i,NS_GlobalConstants); %stimulus information: times at which patterns are played
        [PDChunkNumber,MovieBegin,nRepeats,repeatPeriod,MovieData] = NS_DecodeMovieDataChunk(ChunkData); % Interprets first 6 values in ChunkData.
    end
   
    [Patterns,patternIdx,status] =ReadPatternDataChunk(fnamePattern,PDChunkNumber,NS_GlobalConstants); 
    NumberOfPatterns             =length(patternIdx);
    if i>30
    figure; 
    % Iterates through the repetitions of the current (i-th) movie for the data to subtract.
    for j=1:nRepeats -1
        cla;
        RawData=int16(originalRawFile.getData(MovieBegin + repeatPeriod*(j-1), repeatPeriod)');
        plot(RawData(1+[1:68 317:380],1:2:end)'); 
        title(['movie ' num2str(i) '; repetition ' num2str(j)]);  
        pause(0.1); 
    end
    end
%     keyboard;
end  
fclose(fid0); 