function [raster_times, raster_marks] = spike_raster(spikes, triggers, base_trigger, varargin)
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
%   raster_size           6          the size of each mark on the raster
%
% 2008-10
%   GDF modification of function made by MG


% SET UP OPTIONAL ARGUMENTS

p = inputParser;
p.addRequired('spikes', @isnumeric)
p.addRequired('triggers', @isnumeric)
p.addRequired('base_trigger', @isnumeric)
p.addParamValue('plot_raster', true, @islogical)
p.addParamValue('raster_size', 6, @isnumeric)
p.parse(spikes, triggers, base_trigger, varargin{:});

% triggers denotes the trigger that begins each cycle of the stimulus
raster_size = p.Results.raster_size;

num_cycles = length(triggers) -1;
duration = triggers(2) - triggers(1);

if p.Results.plot_raster
    figure
    clf
    hold on
end

for cycle = 1:num_cycles
    begin_window = triggers(cycle) - base_trigger;
    end_window = triggers(cycle+1) - base_trigger;
    temp_indices = find(spikes >= begin_window & spikes <= end_window);
    temp_raster_times = spikes(temp_indices);
    temp_raster_times = temp_raster_times - begin_window;
    temp_raster_marks = zeros(1,length(temp_raster_times)) + cycle;
    
    raster_times{cycle} = temp_raster_times;
    raster_marks{cycle} = temp_raster_marks;
    
    if p.Results.plot_raster

        temp_rt = raster_times{cycle}';
        temp_rm = raster_marks{cycle};
        for spike = 1:length(temp_raster_times)
            plot([temp_rt(spike), temp_rt(spike)], [(temp_rm(spike)-0.25) (temp_rm(spike)+0.25)], 'k.', 'MarkerSize', raster_size)
        end
        axis([0 duration -1 (num_cycles + 1) ])
    end
end





