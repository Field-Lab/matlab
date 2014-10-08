function elecResp = createElecRespStruct(elecRespInfo)


%% parameters

interval = 50; %length of time before stimulus is updated, in microseconds

elecResp.names.experiment =      elecRespInfo.experimentName; % e.g. '2008-08-28-4'
elecResp.names.data_path =       elecRespInfo.dataPath;
elecResp.names.artifact_path =   elecRespInfo.artifactPath;
elecResp.names.rrs_ei_path =     [elecRespInfo.analysisPath elecRespInfo.analysisBaseName '.ei'];
elecResp.names.rrs_params_path = [elecRespInfo.analysisPath elecRespInfo.analysisBaseName '.params'];
elecResp.names.rrs_sta_path =    '';  % place-holder for future use
elecResp.names.rrs_short_name =  elecRespInfo.shortName;
elecResp.names.savePath =        elecRespInfo.savePath; %location where figures, etc. will be saved

if isfield(elecRespInfo, 'isSpatioTempProbe')
    isSpatioTempProbe = elecRespInfo.isSpatioTempProbe;
else
    isSpatioTempProbe = false;
end


if ~exist(elecResp.names.data_path, 'file') || ~exist(elecResp.names.rrs_ei_path, 'file')...
        || ~exist(elecResp.names.savePath, 'file')
    warnH = warndlg('One of the required paths is not valid');
    uiwait(warnH)
    
    elecResp = [];
    return
end

elecResp.stimInfo.patternNo = elecRespInfo.patternNo;
if ~isempty(elecRespInfo.pElec)
    elecResp.stimInfo.pElec = elecRespInfo.pElec;
else
    elecResp.stimInfo.pElec = 0;
end



if elecRespInfo.autoMovie
    
    movieNos = [];
    if isnumeric(elecResp.stimInfo.patternNo)
        patternNoString = ['p' num2str(elecResp.stimInfo.patternNo)];
    else
        patternNoString = ['p' elecResp.stimInfo.patternNo];
    end
    if exist([elecResp.names.data_path filesep patternNoString], 'file')
        files = dir([elecResp.names.data_path filesep patternNoString]);
    elseif exist(elecResp.names.data_path, 'file')
        files = dir(elecResp.names.data_path);
    else
        error('invalid data path')
    end
    
    for i = 1:length(files)
        if strfind(files(i).name, patternNoString) == 1
            mIndices = strfind(files(i).name, 'm');
            movieNos = [movieNos str2double(files(i).name(mIndices(end)+1:end))]; %#ok<AGROW>
        end
    end
    
    if ~isfield(elecRespInfo, 'autoType') || strcmp(elecRespInfo.autoType, 'all')
    elseif strcmp(elecRespInfo.autoType, 'odds')
        movieNos = movieNos(mod(movieNos,2) == 1);
    elseif strcmp(elecRespInfo.autoType, 'evens')
        movieNos = movieNos(mod(movieNos,2) == 0);
    else
        error('elecRespInfo.autoType is an unexpected string - error in code')
    end
    
    if isempty(movieNos)
        warndlg(['no movies found for pattern ' patternNoString '; skipping to next pattern'])
        elecResp = [];
        return
    end
    
    elecResp.stimInfo.movieNos = sort(movieNos);

    
    if ~isempty(elecResp.names.artifact_path) && exist(elecResp.names.artifact_path, 'file')
        artMovieNos = [];
        if exist([elecResp.names.artifact_path filesep patternNoString], 'file')
            artFiles = dir([elecResp.names.artifact_path filesep patternNoString]);
        elseif exist(elecResp.names.artifact_path, 'file')
            artFiles = dir(elecResp.names.artifact_path);
        else
            error('invalid artifact path')
        end
        for i = 1:length(artFiles)
            if strfind(artFiles(i).name, patternNoString) == 1
                mIndices = strfind(artFiles(i).name, 'm');
                artMovieNos = [artMovieNos str2double(artFiles(i).name(mIndices(end)+1:end))]; %#ok<AGROW>
            end
        end
        
        if ~isfield(elecRespInfo, 'autoType') || strcmp(elecRespInfo.autoType, 'all')
        elseif strcmp(elecRespInfo.autoType, 'odds')
            artMovieNos = artMovieNos(mod(artMovieNos,2) == 1);
        elseif strcmp(elecRespInfo.autoType, 'evens')
            artMovieNos = artMovieNos(mod(artMovieNos,2) == 0);
        else
            error('elecRespInfo.autoType is an unexpected string - error in code')
        end
        
        if isempty(artMovieNos)
            error('no artifact movies found for specified pattern')
        end
        
        elecResp.stimInfo.artMovieNos = sort(artMovieNos);
        
        % removes any movie numbers or artifact movie numbers that are not
        % in both sets of numbers
        sharedMovieNos = intersect(elecResp.stimInfo.movieNos, elecResp.stimInfo.artMovieNos);
        elecResp.stimInfo.movieNos = sharedMovieNos;
        elecResp.stimInfo.artMovieNos = sharedMovieNos;
    else
        elecResp.stimInfo.artMovieNos = [];
    end
    
else
    firstMovie = elecRespInfo.movieFirst;
    lastMovie =  elecRespInfo.movieLast;
    movieInt =   elecRespInfo.movieInt;

    elecResp.stimInfo.movieNos = firstMovie:movieInt:lastMovie;
    
    firstArtMovie = elecRespInfo.artMovieFirst;
    lastArtMovie =  elecRespInfo.artMovieLast;
    artMovieInt =   elecRespInfo.artMovieInt;
    
    if all([~isnan(firstArtMovie) ~isnan(lastArtMovie) ~isnan(artMovieInt)])
        elecResp.stimInfo.artMovieNos = firstArtMovie:artMovieInt:lastArtMovie;
    else
        elecResp.stimInfo.artMovieNos = [];
    end
    
end

nMovies = length(elecResp.stimInfo.movieNos);



elecResp.cells.main = elecRespInfo.mainNeuron;

%includes cells marked as active but NOT main cell
elecResp.cells.active = cell(nMovies,1);

for i = 1:nMovies
    elecResp.cells.active{i} = [];
end

%adds cells to 'active' if specified
for i = 1:4
    x = elecRespInfo.activeNeurons{i};
    if ~isnan(x)
        for j = 1:nMovies
            elecResp.cells.active{j} = [elecResp.cells.active{j} x];
        end
    end
end


elecResp.cells.recElec = elecRespInfo.mainElec;

%.goodElecs includes recElec
elecResp.cells.goodElecs = [];
for i = 1:5
    x = elecRespInfo.otherElecs{i};
    if ~isnan(x)
        elecResp.cells.goodElecs = [elecResp.cells.goodElecs x];
    end
end
elecResp.cells.goodElecs = [elecResp.cells.recElec elecResp.cells.goodElecs];


elecResp.details.sample_rate = elecRespInfo.sampleRate;



% stimulus information determined from given information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elecResp.stimInfo.pulseVectors = cell(nMovies, 1);
elecResp.stimInfo.stimAmps = zeros(nMovies, 1);



%interval = 50; %length of time before stimulus is updated, in microseconds
for i = 1:nMovies
    try
        [amps elecResp.stimInfo.electrodes stimAmpVectors] = getStimAmps(elecResp.names.data_path,...
            elecResp.stimInfo.patternNo, elecResp.stimInfo.movieNos(i));
    catch
       disp('stimulus information couldn''t be retrieved')
       elecResp = [];
       return
    end
    
    %if a primary electrode is specified, stimAmps should reflect amplitude on this electrode
    if any(elecResp.stimInfo.electrodes == elecResp.stimInfo.pElec)
        if ~isSpatioTempProbe
            elecResp.stimInfo.stimAmps(i) = amps(elecResp.stimInfo.electrodes == elecResp.stimInfo.pElec);
        else %stimAmps should reflect "stim pulse" rather than prepulse
            primStimAmpVect = stimAmpVectors(elecResp.stimInfo.electrodes == elecResp.stimInfo.pElec,:);
            stimPulseInd = find(diff(primStimAmpVect~=0)==1, 1, 'last')+1; %index within primStimAmpVect where stim (final) pulse starts
            elecResp.stimInfo.stimAmps(i) = max(abs(primStimAmpVect(stimPulseInd:end)));
        end
    else
        maxIndex = find(abs(amps) == max(abs(amps)),1);
        elecResp.stimInfo.stimAmps(i) = amps(maxIndex);
    end
    
    nElec = length(elecResp.stimInfo.electrodes);
    patternLength = size(stimAmpVectors, 2);
    elecResp.stimInfo.pulseVectors{i} = zeros(nElec, 2, (patternLength-1)*2);
    for j = 1:nElec
        for k = 1:patternLength-1
            elecResp.stimInfo.pulseVectors{i}(j, 1, 2*(k-1)+2) = k*interval;
            elecResp.stimInfo.pulseVectors{i}(j, 1, 2*(k-1)+3) = k*interval;
            elecResp.stimInfo.pulseVectors{i}(j, 2, 2*(k-1)+1) = stimAmpVectors(j,k);
            elecResp.stimInfo.pulseVectors{i}(j, 2, 2*(k-1)+2) = stimAmpVectors(j,k);
        end
        elecResp.stimInfo.pulseVectors{i}(j, 2, 2*(patternLength-1)+1) = stimAmpVectors(j,patternLength);
    end
end

elecResp.stimInfo.nPulses = zeros(1, nMovies);
% determines number of pulses in each movie
for i = 1:nMovies
    dataTraces=NS_ReadPreprocessedData(elecResp.names.data_path, '', 0, elecResp.stimInfo.patternNo,...
        elecResp.stimInfo.movieNos(i));
    elecResp.stimInfo.nPulses(i) = size(dataTraces, 1);
end

% find all neurons that have signal >10 DAQ on any of the center or surrounding electrodes
if isfield(elecRespInfo,'stimSystem')
    switch elecRespInfo.stimSystem
        case '512'
            electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(501);
        case '61'
            electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
    end
else
electrodeMap=edu.ucsc.neurobiology.vision.electrodemap.ElectrodeMapFactory.getElectrodeMap(1);
end
channelsToUse = electrodeMap.getAdjacentsTo(elecResp.cells.recElec, 1)';

for i = 1:length(elecResp.cells.goodElecs)
    if ~any(channelsToUse == elecResp.cells.goodElecs(i))
        channelsToUse = [channelsToUse elecResp.cells.goodElecs(i)]; %#ok<AGROW>
    end
end
elecResp.cells.all = findNeuronsAboveThresh(elecResp.names.rrs_ei_path,...
    elecResp.names.rrs_params_path, channelsToUse, 10);


% check if specified neurons exist in params file
datarun.names.rrs_params_path = elecResp.names.rrs_params_path;
datarun = load_params(datarun, 'verbose', false, 'cell_type_depth', 2, 'sync_cell_ids', true);

if ~any(datarun.cell_ids == elecResp.cells.main) && ~elecRespInfo.externalEi
    warnH = warndlg('Neuron ID specified doesn''t exist in the params file');
    uiwait(warnH)
    elecResp = [];
    return
end

% adds main and active cells to 'all' list if not in it already, and checks to make sure each exists
% in the params file
for i = 1:length(elecResp.cells.active{1})
    if ~any(datarun.cell_ids == elecResp.cells.active{1}(i))
        warnH = warndlg('Neuron ID specified doesn''t exist in the params file');
        uiwait(warnH)
        elecResp = [];
        return
    elseif ~any(elecResp.cells.all == elecResp.cells.active{1}(i))
        elecResp.cells.all = [elecResp.cells.all elecResp.cells.active{1}(i)];
    end
end
if ~any(elecResp.cells.all == elecResp.cells.main)
    elecResp.cells.all = [elecResp.cells.all elecResp.cells.main];
end


% retrieves ei data for cells in elecResp.cells.all
eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(elecResp.names.rrs_ei_path);
if isfield(elecRespInfo,'stimSystem')
    switch elecRespInfo.stimSystem
        case '61'
            numChans = 64;
        case '512'
            numChans = 512; 
    end
else
    numChans = 64; %default for the 61-electrode system
end
if ~elecRespInfo.externalEi
    elecResp.cells.mainEI = eiFile.getImage(elecResp.cells.main);
    elecResp.cells.mainEI = reshape(elecResp.cells.mainEI(1, 2:end, :), numChans, []);
else
    load([elecResp.names.data_path filesep 'eiFile' num2str(elecResp.cells.main)])
    elecResp.cells.mainEI = ei;
end

elecResp.cells.allEIs = cell(length(elecResp.cells.all),1);
for i = 1:length(elecResp.cells.all)
    if elecResp.cells.all(i) ~= elecResp.cells.main
        elecResp.cells.allEIs{i} = eiFile.getImage(elecResp.cells.all(i));
        elecResp.cells.allEIs{i} = reshape(elecResp.cells.allEIs{i}(1, 2:end, :), numChans, []);
    else
        elecResp.cells.allEIs{i} = elecResp.cells.mainEI;
    end
end

eiFile.close()

%since 'active' neurons for all movies should be the same, just check the first movie
for i = 1:length(elecResp.cells.active{1})
    if ~any(datarun.cell_ids == elecResp.cells.active{1}(i))
        warnH = warndlg('One or more of the neuron IDs specified doesn''t exist in the params file');
        uiwait(warnH)
        elecResp = [];
        return
    end
end


% don't need to set (will be set when analyzing) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elecResp.analysis.type =        cell(nMovies,1);
elecResp.analysis.estArtifact = cell(nMovies,1);
elecResp.analysis.latencies =   cell(nMovies,1);
elecResp.analysis.otherLatencies = cell(nMovies,1);
elecResp.analysis.finalized =   zeros(nMovies,1);

elecResp.analysis.successRates = zeros(nMovies,1);
elecResp.analysis.erfParams =   zeros(2,1);
elecResp.analysis.threshold =   0;
elecResp.analysis.erfCurrent =  0; % flag to keep track of whether erfParams and threshold since changing any values in successRates

elecResp.analysis.details.residCalcWindow =  cell(nMovies,1);
elecResp.analysis.details.tempOffsetWindow = cell(nMovies,1);
elecResp.analysis.details.tempOffsetStep =   cell(nMovies,1);
elecResp.analysis.details.linkCalcWindow =   cell(nMovies,1);
elecResp.analysis.details.linkThresh =       cell(nMovies,1);
elecResp.analysis.details.nBranches =        cell(nMovies,1);

elecResp.analysis.details.analysisFlags =    cell(nMovies,1);


for i = 1:nMovies
    elecResp.analysis.details.analysisFlags{i} = zeros(elecResp.stimInfo.nPulses(i), length(elecResp.cells.active{i}) + 1);
end

