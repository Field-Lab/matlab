function subtractArtRawData(dataFolderIn,artDataFolderIn, dataFolderOut,labviewOutputFileFolder, visionWritePath)

% This function removes the electrical artifact in raw .bin files by
% subtracting TTX runs. The raw data and the artifact data must have used
% the same stimulus files.
%
% subtractArtRawData(...) will read the .bin files in dataFolderIn and
% artDataFolderIn, check that the same stimulus patterns were applied to
% both data sets, then average the artifact recording and subtract if from
% the original data. Modified bin files will be created in dataFolderOut
% Lauren Grosberg 9/2014, using WriteDataFile.jar by Georges Goetz

%% Override function inputs for testing
% labviewOutputFileFolder = '/Volumes/Data/2014-08-13-0/';
% dataFolderOut= '/Volumes/Analysis/2014-08-13-0/data003-modified';
% dataFolderIn= '/Volumes/Data/2014-08-13-0/data003';
% artDataFolderIn= '/Volumes/Data/2014-08-13-0/data010';
% visionWritePath = '/Users/grosberg/code/write data GG/WriteDataFile.jar';
if nargin==4
    visionWritePath = '/Users/grosberg/code/write data GG/WriteDataFile.jar';
end

javaaddpath(visionWritePath);

% Link to the relevant files
originalRawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(dataFolderIn);
refHeader = originalRawFile.getHeader();
ttxRawFile = edu.ucsc.neurobiology.vision.io.RawDataFile(artDataFolderIn);
if ~exist(dataFolderOut,'dir')
    mkdir(dataFolderOut);
end
newFile = edu.ucsc.neurobiology.vision.io.ModifyRawDataFile(dataFolderOut, refHeader);

% Getting file info
totalSamples = refHeader.getNumberOfSamples;
Channels = 1:(refHeader.getNumberOfElectrodes-1);
% NS_GlobalConstants = NS_GenerateGlobalConstants(length(Channels));

% Generates the full filenames for movie and pattern output files.
if strcmp(dataFolderIn(end),filesep)
    ix = find(dataFolderIn == filesep,2,'last');
    datarun = dataFolderIn(ix(1)+5:ix(2)-1);
else
    ix = find(dataFolderIn == filesep,1,'last');
    datarun = dataFolderIn(ix(1)+5:end);
end
if strcmp(artDataFolderIn(end),filesep)
    ix = find(artDataFolderIn == filesep,2,'last');
    ttxrun = artDataFolderIn(ix(1)+5:ix(2)-1);
else
    ix = find(artDataFolderIn == filesep,1,'last');
    ttxrun = artDataFolderIn(ix(1)+5:end);
end

fnameMovie      = [labviewOutputFileFolder 'movie' datarun];
% fnamePattern    = [labviewOutputFileFolder 'pattern' datarun];
fnameMovieTtx   = [labviewOutputFileFolder 'movie' ttxrun];
% fnamePatternTtx = [labviewOutputFileFolder 'pattern' ttxrun];

fid0 = fopen(fnameMovie,'r','b');       % Open the movie file.
header = readMHchunk(fid0);             % Read through the header.
nMovies = header.number_of_movies;      % # movie chunks x # amplitudes

fidTtx = fopen(fnameMovieTtx,'r','b');  % Open the movie file.
header = readMHchunk(fidTtx);           % Read through the header.
nMoviesTtx = header.number_of_movies;   % # movie chunks x # amplitudes

if nMovies ~= nMoviesTtx
    err = MException('MATLAB:ArtSubtract:ArtAppDiff', ...
        '# movie chunks X # amplitudes of TTX run does not match original data run');
    throw(err);
end

%%
nSampsPerBinFile = 2400000;
nSampsModified = 0;
tic

for n=1:nMovies-1 % the last movie is expected to be empty
    % Read in the movie parameters for both runs.
    [movieFileChunkNum,movieBegin,nReps,repPeriod,MovieData] = getMovieData(fnameMovie, n);
    [~,movieBeginT,nRepsT,repPeriodT,MovieDataT] = getMovieData(fnameMovieTtx, n);
    
    if ~all(MovieData == MovieDataT)
        err = MException('MATLAB:ArtSubtract:ArtAppDiff', ...
            'TTX pattern application does not match original data');
        throw(err);
    end
    
    %%% Calculate the artifact for the nth movie by averaging repetitions
    try
        ttxData = int16(ttxRawFile.getData(movieBeginT, repPeriod*nRepsT)');
        ttxTraces = reshape(ttxData,size(ttxData,1),[],nRepsT);
        ttxTraces = permute(ttxTraces,[3 1 2]);
    catch
        ttxTraces=int16(zeros(nRepsT,length(Channels)+1, repPeriodT));
        for j=1:nRepsT
            ttxData=int16(ttxRawFile.getData(movieBegin + repPeriod*(j-1), repPeriod)');
            ttxTraces(j,:,:) = ttxData; %RawData(Channels+1,:);
        end
        
    end
    clear ttxData;
    art=int16(squeeze(mean(ttxTraces)))';
    
    % Writing the modified the data to the new file by chunks = repetition
    % period (in samples)
    startSample = movieBegin;
    nSamplesToMod = repPeriod*nReps;
    while nSamplesToMod > 0
        data = originalRawFile.getData(startSample, repPeriod);
 
        modifiedData = int16(zeros(repPeriod, length(Channels)+1)); 
        modifiedData(:,1) = data(:,1); %Leave TTL pulses alone
        modifiedData(:, 2:end) = data(:, Channels+1) - art(:, Channels+1);
        nSampsModified = nSampsModified + min(repPeriod,nSamplesToMod);
        if nSampsModified > nSampsPerBinFile
            nSampsModified = 0;
            newFile.addFile();
        end
        newFile.appendDataToLastFile(modifiedData);
        nSamplesToMod = nSamplesToMod - repPeriod;
        startSample = startSample + repPeriod;
    end
    
     fprintf('finished processing %d of %d movies in %0.1f s. estimated time remaining: %0.1f min\n',n,nMovies-1,toc, toc/n*(nMovies-1-n)/60);
     
end
toc
fprintf('Finished \n');
% Closing the files
originalRawFile.close();
newFile.close();