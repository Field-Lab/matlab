function neuronsAboveThresh = findNeuronsAboveThresh(eiPath, paramsPath, channels, thresh, neurons)
% usage: neuronsAboveThresh = findNeuronsAboveThresh(eiPath, paramsPath, channels, thresh, [neurons])
%
% gets a list of all neurons in params file that have an ei with signal (absolute value) above the
% specified threshold on at least one of the specified channels
%
% arguments
%   eiPath: full path to .ei file
%   paramsPath: full path to .params file
%   channels: electrodes that you want to check for > thresh DAQ signal on
%   thresh: threshold DAQ value of an ei, on one of the chosen channels, that's required for the
%   neuron to be included in the returned list
%   neurons (optional): list of neurons to check (default = all neurons from associated
%   params file that are not classified as duplicates or contaminated)
%
% returns
%   neuronsAboveThresh: vector of neuron IDs that have signal above threshold on at least one of the
%   specified electrodes
%
% author: Lauren Hruby (lhruby@salk.edu)
% last updated: 2009-06-17
%


datarun.names.rrs_params_path = paramsPath;
params.req_order = {'dup','contam','contam '};
params.avoid = '';

datarun = load_params(datarun, 'verbose', false, 'cell_type_depth', 2, 'sync_cell_ids', true, 'order_cell_types', false);
datarun.cell_types = order_cell_types(datarun.cell_types, params);

if exist('neurons', 'var')
    for i = length(neurons):-1:1
        %remove any cell ids that don't exist in params file
        if ~any(datarun.cell_ids == neurons(i)) %if neuron in chosen list is not in .params file
            warnH = warndlg('One or more of the neuron IDs specified doesn''t exist in the params file');
            uiwait(warnH)
            neurons(i) = [];
        %remove any contaminated or duplicate cell ids from list
        elseif any(datarun.cell_types{1}.cell_ids == neurons(i)) ||...
                any(datarun.cell_types{2}.cell_ids == neurons(i)) ||...
                any(datarun.cell_types{3}.cell_ids == neurons(i))
            neurons(i) = [];
        end
    end
else %use all neurons in .params file
    neurons = [];
    for i = 1:length(datarun.cell_ids)
        if ~any(datarun.cell_types{1}.cell_ids == datarun.cell_ids(i)) &&...
                ~any(datarun.cell_types{2}.cell_ids == datarun.cell_ids(i)) &&...
                ~any(datarun.cell_types{3}.cell_ids == datarun.cell_ids(i))
            neurons = [neurons datarun.cell_ids(i)]; %#ok<AGROW>
        end
    end
end


nNeurons = length(neurons);
neuronMaxs = zeros(nNeurons, 1);

eiFile = edu.ucsc.neurobiology.vision.io.PhysiologicalImagingFile(eiPath);
for i = 1:nNeurons
    eiData = eiFile.getImage(neurons(i));
    eiData = squeeze(eiData(1, 2:end, :));
    
    neuronMaxs(i) = max(max(abs(eiData(channels, :))));
end
eiFile.close()

neuronsAboveThresh = neurons(neuronMaxs >= thresh);