function datarun = get_interspikeinterval(datarun, cell_spec, varargin)
% MY_FUNCTION    compute the inter-spike-interval for all cells specified in
%                cell_spec
%
% usage:  datarun = get_interspikeinterval(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct 
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field
%                       interspikeinterval.probabilities{}
%                       interspikeinterval.bins{}    
%
% optional parameters, their default values, and what they specify:
%
%
% bin_size       0.001              bin size for autocorrelogram (in
%                                   seconds)
% duration        0.1               duration to calculate autocorrelation
%                                   (in seconds)       
%
%
%
% 2011-02  GDF
% 2014-04  sravi modified function name and function call to interspikeintervals
%          sravi added ISI calculation per trial for pulses and drifting gratings  


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('bin_size', 0.001, @isnumeric);
p.addParamValue('duration', 0.1, @isnumeric);
p.addParamValue('pulses', false);
p.addParamValue('dg', false);
p.addParamValue('sp', 64, @isnumeric);
p.addParamValue('tp', 256, @isnumeric);
p.addParamValue('trial_begin_times', [], @isnumeric);
p.addParamValue('stop', [], @isnumeric);


% resolve user input and default values
p.parse(varargin{:});

% BODY OF FUNCTION

% get number of cells and indices
cell_indices = get_cell_indices(datarun, cell_spec);
num_cells = length(cell_indices);

% initialize interspikeinterval field if it does not already exist
if ~isfield(datarun, 'interspikeinterval')
    datarun.interspikeinterval = cell(length(datarun.cell_ids),1);
end

% get interspikeinterval and bins for cells of interest
for cc = 1:num_cells
    
    spike_times = datarun.spikes{cell_indices(cc)};
    
    if(p.Results.pulses)
        spikes_by_trials = get_raster(spike_times, p.Results.trial_begin_times, 'plot', false);
        
        for tr = 1:length(spikes_by_trials)
            [prob, bi] = interspikeinterval(spikes_by_trials{tr,1}, p.Results.bin_size, p.Results.duration, mean(diff(p.Results.trial_begin_times)));
            probabil(tr,:) = prob;
        end
        
        probabilities = sum(probabil)./(length(spikes_by_trials));
        bins = bi;
    
    elseif(p.Results.dg)
        for i = 1:length(datarun.stimulus.params.DIRECTION)
            trnum = intersect(intersect(find([datarun.stimulus.combinations.SPATIAL_PERIOD] == p.Results.sp),find([datarun.stimulus.combinations.TEMPORAL_PERIOD] == p.Results.tp)), find([datarun.stimulus.combinations.DIRECTION] == datarun.stimulus.params.DIRECTION(1,i)));
            trigpre = ismember(datarun.stimulus.trial_list,trnum);
            spikes_by_trials = get_raster(spike_times, datarun.stimulus.triggers(trigpre), 'stop', mean(diff(datarun.stimulus.triggers)),'plot', false);
            
            for tr = 1:length(spikes_by_trials)
                [prob, bi] = interspikeinterval(spikes_by_trials{tr,1}, p.Results.bin_size, p.Results.duration,  mean(diff(datarun.stimulus.triggers)));
                probabil(tr,:) = prob;
            end
            
            probabilit(i,:) = sum(probabil)./(length(spikes_by_trials));
        end
        
        probabilities = sum(probabilit)./(length(datarun.stimulus.params.DIRECTION));
        bins = bi;         
        
    else
        [probabilities, bins] = interspikeinterval(spike_times, p.Results.bin_size, p.Results.duration, datarun.duration);
    end

    % store information in datarun
    datarun.interspikeinterval{cell_indices(cc)}.probabilities = probabilities;
    datarun.interspikeinterval{cell_indices(cc)}.bins = bins;
end

end

    