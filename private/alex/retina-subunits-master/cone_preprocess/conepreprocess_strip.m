function stripped = conepreprocess_strip(datarun)
stripped = datarun;

if isfield(stripped, 'spikes')
    stripped = rmfield(datarun, 'spikes');
end

if isfield(stripped.stas, 'stas')
    stripped.stas = rmfield(stripped.stas, 'stas');
end

if isfield(stripped, 'spike_rate')
    stripped.spike_rate = uint8(stripped.spike_rate);
end