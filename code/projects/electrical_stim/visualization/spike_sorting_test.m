% Script to analyze data using Gonzalo's spike sorting algorithm.
codebase_path = matlab_code_path;
analysisPath = uigetdir('/Volumes/Analysis/', 'Choose the electrical stim data folder (001,002, etc) that contains data organized by pattern');
[eiFileName, eiPath] = uigetfile('/Volumes/Analysis/*.ei','Choose ei file to use as a template');
pathToEi = [eiPath eiFileName];

patternNos = input('Enter pattern(s) to create sorted output file(s): ');
% neuronIds = input('Enter neuron ids to create sorted output file(s): ');
% analysisPath = '/Volumes/Analysis/2015-05-27-0/data001/'; % Directory with the electrical spikes sorted into -autosort/';

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
sortedPath = [analysisPath '-autosort'];
if ~isdir(sortedPath)
    mkdir(sortedPath);
end

for n = 1:length(neuronIds)
    % Save output files in the created directory
    for p = 1:length(patternNos)
        patternNo = patternNos(p);
        elecRespAuto = SpikeSortingCompact(analysisPath,patternNo,neuronIds(n),pathToEi,...
            'recElecs',recElecs,'findAxon',findAxon,'findAxon',findAxon,...
            'cleanData',cleanData,'Trange',Trange,'degPolRule',degPolRule);
        fname = fullfile(sortedPath,['elecRespAuto_n' ...
            num2str(elecRespAuto.neuronInfo.neuronIds) '_p' ...  % Should be neuronIds(n)?
            num2str(elecRespAuto.stimInfo.patternNo) '.mat']);
        save(fname,'elecRespAuto');
        disp(['done analyzing ' fname]);
        %     save([elecResp.names.data_path filesep elecRespName], 'elecResp')
        %         disp(['done analyzing movie ' num2str(elecResp.stimInfo.movieNos(j)) ', pattern ' num2str(patternNos(i))])
    end
    
end



%translate templates if by some reason there are weird template offsets.
%templates = translateTemplate(templates,Translate(p),1,1); % The minimum of each
% template should always align with sample point 10