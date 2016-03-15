function datarun = convert_datarun_times_to_samples(datarun)
%CONVERT_DATARUN_TIMES_TO_SAMPLES Converts all the datarun times to
%samples. This function assumes the datarun times were previously specified
%in seconds.

if ~isfield(datarun, 'sampling_rate')
    error('Sampling rate has not yet been loaded')
end

if isfield(datarun, 'triggers')
    datarun.triggers = round(datarun.triggers * datarun.sampling_rate);
end

if isfield(datarun, 'duration')
    datarun.duration = round(datarun.duration * datarun.sampling_rate);
end

if isfield(datarun, 'spikes')
    for k = 1:length(datarun.spikes)
        datarun.spikes{k} = round(datarun.spikes{k} * datarun.sampling_rate);
    end
end

end % convert_datarun_times_to_samples