function combineERFMovies(outputFolderStatistics,allSourceFolders)
% This function combines the eRFs movies computed for different pixel 
% positions into a single eRF movie structure, and plots the result.
%
% Parameters:
%   - outputFolderStatistics: statistics will be saved here
%   - allSourceFolders: a cell with all source statistics folders
%
% Returns:
%   []
%

if nargin<2
    error('Matlab:PlotERFs:combineAndPlotERFs','Please provide at least one source statistics folder');
end
% allSourceFolders = allSourceFolders{1};
nSourceFolders = length(allSourceFolders);


%% Formatting the folders 

if outputFolderStatistics(end:end)~=filesep
    outputFolderStatistics = [outputFolderStatistics filesep];
end
for kk=1:nSourceFolders
    if allSourceFolders{kk}(end:end)~=filesep
        allSourceFolders{kk} = [allSourceFolders{kk} filesep];
    end
end

% If output folders don't exist, make them here
if ~exist(outputFolderStatistics,'dir')
    mkdir(outputFolderStatistics);
end

%% Building the neuronlist

neuronList = [];
for kk=1:nSourceFolders
    neuronList = unique(union(neuronList, readNeuronNames(allSourceFolders{kk})));
end
nNeurons = length(neuronList);

%% Combine the neurons eRFs

for kk=1:nNeurons
    neuronTdERFdata = combineNeuronERFs(neuronList(kk), allSourceFolders);

    % Save the eRF
    save(fullfile(outputFolderStatistics, sprintf('neuron%d_td_erf',neuronList(kk))),'neuronTdERFdata');
    
end

end % combineERFMovies

function neuronList = readNeuronNames(statsFolder)
% This function returns a list of all the neuron IDs found in a statistics
% folder.

contentsActivationStatsFolder = dir(statsFolder);
neuronList = [];

for kk=1:length(contentsActivationStatsFolder)
    if strfind(contentsActivationStatsFolder(kk).name,'.mat')
        neuronList(end+1) = str2double(contentsActivationStatsFolder(kk).name(7:end-11));
    end
end

end

function neuronTdERFdatacombined = combineNeuronERFs(neuronID, allSourceFolders)
% This function combines all the erf data into a single structure.

% Get an existing structure from one of the source folders and load it.
for kk=1:length(allSourceFolders)
    try load(fullfile(allSourceFolders{kk},sprintf('neuron%d_td_erf',neuronID)));
        break;
    catch %#ok<CTCH>
    end
end
neuronTdERFdatacombined = neuronTdERFdata;

% Now that we've got this stucture, fill it in
for kk=1:length(allSourceFolders)
    try load(fullfile(allSourceFolders{kk},sprintf('neuron%d_td_erf',neuronID)));
        neuronTdERFdatacombined.allPSTHs(neuronTdERFdata.pixelIndex,:) = ...
            neuronTdERFdata.allPSTHs(neuronTdERFdata.pixelIndex,:);
        neuronTdERFdatacombined.pixelIndex = union(neuronTdERFdatacombined.pixelIndex,...
            neuronTdERFdata.pixelIndex);
    catch %#ok<CTCH>
    end
end

end