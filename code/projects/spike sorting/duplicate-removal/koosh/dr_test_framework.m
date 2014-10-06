% dr_test_framework
%
% this wrapper script reads in the datarun specified and uses matlab to do
% neuron cleaning using the neurons file and ei files specified below. note
% that the function will reversibly modify the input neuron file by
% flagging neurons as "removed". consequently, you should make a copy of
% the raw neurons file and possibly add a suffix to it to differentiate it
% from the original raw neurons file.
%
clear all; close all;

% get the sampling rate
samplingRate = edu.ucsc.neurobiology.vision.util.VisionParams.samplesPerMillisecond * 1000;

% specify the dataset analysis folder
analysisFolder = '/Data/Machado/2009-04-13-0/data000/test-dr/test-dr';

% initialize datarun
datarun = load_data(analysisFolder);

% update the path to the neurons-raw file:
if 1
    % name of the neurons file (that will be modified)
    datarun.names.rrs_neurons_path = [datarun.names.rrs_neurons_path '-copy'];
    % name of input ei file (computed on ALL cells in neurons-raw file)
    datarun.names.rrs_ei_path = [datarun.names.rrs_ei_path];
end

% load neurons
datarun = load_neurons(datarun);

% get cell numbers
indices = get_cell_indices(datarun, datarun.cell_ids);

% get spike trains for all cells that we want to add to the neurons file
processed_cells = dr_process_duplicates(datarun, samplingRate);

% update the neurons file:
if 1
    % open neurons file
    nf = edu.ucsc.neurobiology.vision.io.NeuronFile(datarun.names.rrs_neurons_path);

    % get rid of all neurons
    kill = [];
    for ii = 1:length(datarun.cell_ids)
        nf.deleteNeuron(datarun.cell_ids(ii));
    end

    % add neurons
    id = 0;
    for ii = 1:length(processed_cells.s)
        il = edu.ucsc.neurobiology.vision.util.IntegerList;
        spikes = processed_cells.s{ii};
        for jj = 1:length(processed_cells.s{ii})
            il.add(spikes(jj));
        end
        
        nf.addNeuron(processed_cells.e(ii), id, il,length(processed_cells.s{ii}))
        id = id + 1;
    end
    
    nf.close;
end


    
