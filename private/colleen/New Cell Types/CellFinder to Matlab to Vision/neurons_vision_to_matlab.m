% neurons_vision_to_matlab
% Colleen Rhoades
% June 2015
% rhoades@stanford.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The purpose of this file is to load an neurons file and write a new 
% neurons files with just a subset of those neurons.

% This is useful for being able to export all neurons from Cell Finder and 
% then cut the file down to just those you are interested in to load in 
% Vision.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

date = '2016-02-17-6';
datarun = '025'; % just last three numbers

% Data that contains all the neurons
referenceRawDataFilePath = ['/Volumes/Data/' date, '/data', datarun,'/'];

% File that contains all the neurons or params (saved from Cell Finder or computed by
% Grind)
% referenceNeuronFilePath = ['/Volumes/Analysis/', date, '/data011/data011_allneurons/data011_allneurons_manual/data011_allneurons_manual.neurons'];
% referenceParamsFilePath = ['/Volumes/Analysis/', date, '/data011/data011_allneurons/data011_allneurons_manual/data011_allneurons_manual.params'];

referenceNeuronFilePath = ['/Volumes/Analysis/', date, '/data', datarun,'-cf/data', datarun,'-cf.neurons'];
referenceParamsFilePath = ['/Volumes/Analysis/', date, '/data', datarun,'-cf/data', datarun,'-cf.params'];


% Location where a new neurons file containing a subset of the neurons will
% be located
newNeuronFilePath = ['/Volumes/Analysis/', date, '/data', datarun,'-cf/edited/data', datarun,'-cf/data', datarun,'-cf.neurons'];

% Set specify_cell to 1 if you want to list cell ids. If you want to get a
% whole class of cells (say all ON Parasols) then set specify_cells to 0
specify_cells = 1;

% Specify which cell class you want (can be anything if specify_cells = 1)
% cell_type = {'OFF parasol'};
% 
% if specify_cell = 0, list cell ids here (these are Vision IDs)
toKeep = [
 1216 1936 2011 3901 4111 4981 5013 5523 5646 6111 6528 6953 7051 7457 7471
    ];
toKeep = sort(toKeep);
%%%%%%%%%%%%%%%%%%%%%%%%%% END INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Data
rawDataFile = edu.ucsc.neurobiology.vision.io.RawDataFile(referenceRawDataFilePath);
oldNeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(referenceNeuronFilePath);
ttlTimes = oldNeuronFile.getTTLTimes();

% Run a function Georges wrote to get the right header
newNeuronFile = initVisionNeuronFile(referenceRawDataFilePath, newNeuronFilePath, ttlTimes);


if specify_cells == 0
    datarun.names.rrs_neurons_path=[referenceNeuronFilePath];
    datarun.names.rrs_params_path=[referenceParamsFilePath];
    opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false, 'load_sta', 0, 'load_sta_params', 0, 'load_all',false);
    datarun=load_data(datarun,opt);
    
    
    cell_type_index= zeros(1,size(cell_type,2));
    for num_cell_types = 1:size(cell_type,2)
        for i = 1:size(datarun.cell_types,2)
            right_cell_type = strcmpi(datarun.cell_types{i}.name, cell_type{num_cell_types}); % case insensitive
            if right_cell_type == 1;
                cell_type_index(num_cell_types) = i;
                break
            end
            cell_type_index(num_cell_types) = 0;% couldn't find the right cell type
        end
        
    end
    cell_specification = datarun.cell_types{cell_type_index}.cell_ids;
    % The cell_specification is all ids from the specified class
    toKeep = [
        cell_specification
        ];
    
end

% Write each selected neuron to the new neurons file.
for i = 1:length(toKeep)
    spikeTimes = oldNeuronFile.getSpikeTimes(toKeep(i));
    electrode = oldNeuronFile.getNeuronIDElectrode(toKeep(i));
    newNeuronFile.addNeuron(electrode, toKeep(i), spikeTimes, length(spikeTimes));
end

% Close the files so you can open them in Vision
newNeuronFile.close();
oldNeuronFile.close();
rawDataFile.close();

