% Java library
% javaaddpath('/Applications/Vision.app/Contents/Resources/Java/Vision.jar');
javaaddpath('/Volumes/Lab/Development/vision7/Vision.app/Contents/Resources/Java/Vision.jar');

% matlab code paths 
addpath(genpath('/home/vision/Dropbox/Lab/Development/matlab-standard/code'));

% load dataset
datarun = load_data('2005-04-26-0/data000-nwpca/data000');
datarun = load_sta(datarun);
datarun = load_params(datarun);
datarun = load_neurons(datarun);
datarun = set_polarities(datarun);
datarun = load_ei(datarun, 'all');

% get the on midget cell ids
onmidget_ids = get_cell_indices(datarun,{'on on midget nc4'});
onmidget_vision_ids = datarun.cell_ids(onmidget_ids);
N = length(onmidget_ids);

% to get spike train for one cell, use datarun.spikes{onmidget_ids(i)}

% get the on midget cell electrodes
onmidget_electrodes = datarun.channels(onmidget_ids);

% get all cells on these electrodes
all_cells_on_electrodes = [];
for i = 1:N
    all_cells_on_electrodes = ...
        [all_cells_on_electrodes; 
         cell_ids_from_electrode(onmidget_electrodes(i))]; %#ok<*AGROW>
end
all_cells_on_electrodes = unique(all_cells_on_electrodes);

% get waveforms for on midget cells
% first index corresponds to cell, second to electrode, third to timepoint
onmidget_waveforms = nan(N,512,81);
for i = 1:N
    onmidget_waveforms(i,:,:) = datarun.ei.eis{onmidget_ids(i)};
end