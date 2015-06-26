clear
date = '2015-04-14-2';
referenceRawDataFilePath = ['/Volumes/Data/' date, '/data008/'];
rawDataFile = edu.ucsc.neurobiology.vision.io.RawDataFile(referenceRawDataFilePath);

referenceNeuronFilePath = ['/Volumes/Analysis/', date, '/data008-cr/data008.neurons'];
oldNeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(referenceNeuronFilePath);
ttlTimes = oldNeuronFile.getTTLTimes();

newNeuronFilePath = ['/Volumes/Analysis/', date, '/data008-cr/data008-cr.neurons'];
newNeuronFile = initVisionNeuronFile(referenceRawDataFilePath, newNeuronFilePath, ttlTimes);

specify_cells = 1;
if specify_cells == 0
referenceParamsFilePath = ['/Volumes/Analysis/', date, '/data002-gdf/data002/data002.params'];
datarun.names.rrs_neurons_path=[referenceNeuronFilePath];
datarun.names.rrs_params_path=[referenceParamsFilePath];
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false, 'load_sta', 0, 'load_sta_params', 0, 'load_all',false);
datarun=load_data(datarun,opt);

cell_type = {'OFF parasol'};
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
toKeep = [
    cell_specification
    ];
else
    
toKeep = [
1561
2302
3486
5391
7156
7157
396
2242
3035
3361
4853
6076
6993
7486
7487
7488
    ];
end

for i = 1:length(toKeep)
    spikeTimes = oldNeuronFile.getSpikeTimes(toKeep(i));
    electrode = oldNeuronFile.getNeuronIDElectrode(toKeep(i));
    newNeuronFile.addNeuron(electrode, toKeep(i), spikeTimes, length(spikeTimes));
end

newNeuronFile.close();
oldNeuronFile.close();
rawDataFile.close();

