function [raster_times, raster_marks] = spike_raster(spikes, triggers, varargin)
% spike_raster      Creates cell-arrays for generating spike rasters for a
%                   periodic stimulus. It can also generate a plot of the
%                   rasters
%
% usage:   [raster_times, raster_marks] = spike_raster(spikes, triggers, plot_raster)
%
% arguments:    spikes - vector of spike times
%               triggers - vector of triggers denoting beginning of each cycle 
%
% outputs:  raster_times - times of spikes sorted by stimulus cycle  
%           raster_marks - denotes trial or cycle of the stimulus
%
% optional fields with defaults
%   
%   plot_raster         true         will generate a plot of the raster
%
% 2008-10
%   GDF modification of function made by MG


% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addRequired('spikes', @isarray)
p.addRequired('triggers', @isarray)
p.addParamValue('plot_raster', true, @islogical)
p.parse(master_datarun, slave_datarun, varargin{:});

% triggers denotes the trigger that begins each cycle of the stimulus

num_cycles = length(triggers) - 1;
duration = round(triggers(2) - triggers(1));

for cycle = 1:num_cycles
    begin_window = triggers(cycle) - triggers(1);
    end_window = triggers(cycle+1) - triggers(1);
    temp_indices = find(spikes >= begin_window & spikes <= end_window);
    temp_raster_times = spikes(temp_indices);
    temp_raster_times = temp_raster_times - (duration *(cycle-1));
    temp_raster_marks = zeros(1,length(temp_raster_times)) + cycle;
    
    raster_times{cycle} = temp_raster_times;
    raster_marks{cycle} = temp_raster_marks;
    
    if p.Results.plot_raster
        plot(raster_times{cycle}, raster_marks{cycle}, 'k.', 'MarkerSize', 6)
        axis([0 duration -1 (num_cycles + 1) ])
        hold on
    end
end





