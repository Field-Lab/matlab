function datarun = get_autocorrelations(datarun, cell_spec, varargin)
% MY_FUNCTION    compute the autocorrelation for all cells specified in
%                cell_spec
%
% usage:  datarun = get_auto_correlations(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct 
%           cell_spec - which cells - cell ids
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    datarun - datarun struct with results stored in field
%                       autocorrelation.probabilities{}
%                       autocorrelation.bins{}    
%
% optional parameters, their default values, and what they specify:
%
% sampling_freq       2kHz               rate at which to sample spike train
%
% duration            0.1                duration to calculate autocorrelation (in seconds)
% 
% pulses              false              whether to calculate ACF trial by trial for full field pulse stimulus
%
% trial_begin_times   []                 vector of Full Field Pulse trial begin times - must specify when pulses = true
%
% stop                []                 time to stop trial, default is mean of trigger differences
%
% dg                  false              whether to calculate ACF trial by trial for a certain spatial and temporal period within the drifting grating stimulus
%
% sp                  64                 spatial Period of Drifting Grating
%
% tp                  256                temporal Period of Drifting Grating
%
% acf_func            0                  whether to use autocorrelation.m to calculate ACFs (more accurate option but slower) - samples all spikes, xcorr can miss few spikes
%
% verbose             true              display wait bar
%
% 2011-02  GDF
% 2014-01  sravi modified autocorrelation function
%          sravi added ACF calculation per trial for pulses and drifting gratings 


% SET UP OPTIONAL ARGUMENTS

p = inputParser;

% specify list of optional parameters
p.addParamValue('sampling_freq', 2000, @isnumeric);
p.addParamValue('duration', 0.1, @isnumeric);
p.addParamValue('pulses', false);
p.addParamValue('trial_begin_times', [], @isnumeric);
p.addParamValue('stop', [], @isnumeric);
p.addParamValue('dg', false);
p.addParamValue('sp', 64, @isnumeric);
p.addParamValue('tp', 256, @isnumeric);
p.addParamValue('acf_func', 0, @isnumeric);
p.addParamValue('verbose', true, @islogical);

% resolve user input and default values
p.parse(varargin{:});
params = p.Results;

% BODY OF FUNCTION

% get number of cells and indices
cell_indices = get_cell_indices(datarun, cell_spec);
num_cells = length(cell_indices);
bin_size = 1/params.sampling_freq;

% initialize autocorrelation field if it does not already exist
if ~isfield(datarun, 'autocorrelation')
    datarun.autocorrelation = cell(length(datarun.cell_ids),1);
end

% show output
if params.verbose
    T = text_waitbar(sprintf('Getting autocorrelation functions...', length(cell_indices)));
    start_time = clock; % note when it started
end


% get autocorrelation and bins for cells of interest
for cc = 1:num_cells
    
    spike_times = datarun.spikes{cell_indices(cc)};
    
    if(params.pulses)
        spikes_by_trials = get_raster(spike_times, params.trial_begin_times, 'plot', false);
        probabil = zeros(length(params.trial_begin_times), (params.duration/bin_size)+1);
        for tr = 1:length(spikes_by_trials)
            
            if(params.acf_func)
                [prob, bi] = autocorrelation(spikes_by_trials{tr,1}, bin_size, params.duration, mean(diff(params.trial_begin_times)));
                probabil(tr,:) = prob;
            else
                sigVec = 0:1:(params.sampling_freq*round(mean(diff(params.trial_begin_times))));
                spikeTimes = round(spikes_by_trials{tr,1}*(params.sampling_freq));
                spikeVec = double(ismember(sigVec, spikeTimes)); %vector of zeroes and ones of spike instances
                [prob bi] = xcorr(spikeVec,params.sampling_freq*params.duration); %lags need to be an integer
                probabil(tr,:) = [0 prob(params.sampling_freq*params.duration+2:end)];
            end                       
        end
        
        probabilities = sum(probabil);
        bins = 0:bin_size:params.duration;
        
        elseif(params.dg)
        
        for i = 1:length(datarun.stimulus.params.DIRECTION)
            trnum = intersect(intersect(find([datarun.stimulus.combinations.SPATIAL_PERIOD] == params.sp),find([datarun.stimulus.combinations.TEMPORAL_PERIOD] == params.tp)), find([datarun.stimulus.combinations.DIRECTION] == datarun.stimulus.params.DIRECTION(1,i)));
            trigpre = ismember(datarun.stimulus.trial_list,trnum);
            if isempty(params.stop)
                params.stop = mean(diff(datarun.stimulus.triggers));
            end
            spikes_by_trials = get_raster(spike_times, datarun.stimulus.triggers(trigpre), 'stop', params.stop, 'plot', false);
            probabil = zeros(sum(trigpre), (params.duration/bin_size)+1);
            probabilit = zeros(length(datarun.stimulus.params.DIRECTION), (params.duration/bin_size)+1);
            for tr = 1:length(spikes_by_trials)
                if(params.acf_func)
                [prob, bi] = autocorrelation(spikes_by_trials{tr,1}, bin_size, params.duration,  mean(diff(datarun.stimulus.triggers)));
                probabil(tr,:) = prob;
                else
                    sigVec = 0:1:(params.sampling_freq*round(mean(diff(datarun.stimulus.triggers))));
                    spikeTimes = round(spikes_by_trials{tr,1}*(params.sampling_freq));
                    spikeVec = double(ismember(sigVec, spikeTimes)); %vector of zeroes and ones of spike instances
                    [prob bi] = xcorr(spikeVec,params.sampling_freq*params.duration); %lags need to be an integer
                    probabil(tr,:) = [0 prob(params.sampling_freq*params.duration+2:end)];
                end
                    
            end
            
            probabilit(i,:) = sum(probabil);
        end
        
        probabilities = sum(probabilit);
        bins = 0:bin_size:params.duration;        
        
       else
                if(params.acf_func)
                [probabilities, bins] = autocorrelation(spike_times, bin_size, params.duration, datarun.duration);
                else
                    sigVec = 0:1:(params.sampling_freq*datarun.duration);
                    spikeTimes = round(spike_times*(params.sampling_freq));
                    spikeVec = double(ismember(sigVec, spikeTimes)); %vector of zeroes and ones of spike instances
                    [prob bi] = xcorr(spikeVec,params.sampling_freq*params.duration); %lags need to be an integer
                    probabilities = [0 prob(params.sampling_freq*params.duration+2:end)];
                    bins = 0:bin_size:params.duration;        
                end
    end 
    

    % store information in datarun
    datarun.autocorrelation{cell_indices(cc)}.probabilities = probabilities;
    datarun.autocorrelation{cell_indices(cc)}.bins = bins;
    if params.verbose
        T = text_waitbar(T, cc/length(cell_indices));
    end
end

end

    