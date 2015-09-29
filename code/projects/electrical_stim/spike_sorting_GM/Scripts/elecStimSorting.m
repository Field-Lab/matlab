function elecStimSorting(pathToAnalysisData,pathToEi,patternNos,neuronIds)
% inputs:     pathToAnalysisData: '/Volumes/Analysis/...'choose the electrical stim data folder (001,002, etc) that contains data organized by pattern
%             pathToEi: a string, '/full/path/to/ei/file.ei'
%             patternNos
%             neuronIds


% codebase_path = matlab_code_path; 

% Set optional arguments. 
recElecs=[];                % Electrodes that will be used for spike sorting. If no value is specified
                            % then electrodes from the electrode with the largest
                            % EI signal will be used for each neuron
findAxon = 0;               % logical, if zero no axon breakpoints will be included (default =0)
cleanData = 1;              % if 1, the first trial of each movie will not be included. (default =1)
collapseTrials = 1;         % if 1, movies corresponding to same stimulus will be collapsed into one (default =1)
Trange = [1 40];            % sample range over which to do spike sorting, specified as a two dimensional vector (default = [1 40])
degPolRule = [1 4 15 50]; 

% Create directory for automatic spike sorting
if strcmp(pathToAnalysisData(end),filesep)
    sortedPath = [pathToAnalysisData(1:end-1) '-autosort'];
else
    sortedPath = [pathToAnalysisData '-autosort'];
end
if ~isdir(sortedPath)
    mkdir(sortedPath);
end

% Save output files in the created directory
for p = 1:length(patternNos)
    patternNo = patternNos(p); 
    elecRespAuto = SpikeSortingCompact(pathToAnalysisData,patternNo,neuronIds,pathToEi,...
        'recElecs',recElecs,'findAxon',findAxon,'findAxon',findAxon,...
        'cleanData',cleanData,'Trange',Trange,'degPolRule',degPolRule);
    fname = fullfile(sortedPath,['elecRespAuto_n' ...
        num2str(elecRespAuto.neuronInfo.neuronIds) '_p' ...
        num2str(elecRespAuto.stimInfo.patternNo) '.mat']); 
    save(fname,'elecRespAuto'); 
    disp(['done analyzing ' fname]); 
end