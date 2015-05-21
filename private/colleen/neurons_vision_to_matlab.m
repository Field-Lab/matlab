clear

referenceRawDataFilePath = '/Volumes/Archive/2010-08-20-1/data001/';
rawDataFile = edu.ucsc.neurobiology.vision.io.RawDataFile(referenceRawDataFilePath);

referenceNeuronFilePath = '/Volumes/Analysis/2010-08-20-1/data001-nwpca-copy/data001-nwpca/data001-nwpca.neurons';
oldNeuronFile = edu.ucsc.neurobiology.vision.io.NeuronFile(referenceNeuronFilePath);
ttlTimes = oldNeuronFile.getTTLTimes();

newNeuronFilePath = '/Volumes/Analysis/2010-08-20-1/data001-cr/data001-cr.neurons';
newNeuronFile = initVisionNeuronFile(referenceRawDataFilePath, newNeuronFilePath, ttlTimes);

referenceParamsFilePath = '/Volumes/Analysis/2010-08-20-1/data001-nwpca-copy/data001-nwpca/data001-nwpca.params';
datarun.names.rrs_neurons_path=[referenceNeuronFilePath];
datarun.names.rrs_params_path=[referenceParamsFilePath];
opt=struct('verbose',1,'load_params',1,'load_neurons',1,'load_obvius_sta_fits',false, 'load_sta', 0, 'load_sta_params', 0, 'load_all',false);
datarun=load_data(datarun,opt);

cell_type = {'OFF large 2'};
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

for i = 1:length(toKeep)
    spikeTimes = oldNeuronFile.getSpikeTimes(toKeep(i));
    electrode = oldNeuronFile.getNeuronIDElectrode(toKeep(i));
    newNeuronFile.addNeuron(electrode, toKeep(i), spikeTimes, length(spikeTimes));
end

newNeuronFile.close();
oldNeuronFile.close();
rawDataFile.close();

