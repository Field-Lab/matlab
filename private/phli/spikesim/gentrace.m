function trace = gentrace(waveforms, sts, varargin)
% GENTRACE
% usage: gentrace(waveforms, sts, opts)
%
% Spike times STS is assumed to be in sorted order.  Waveforms is cell
% array of waveforms, which can be vectors or matrices.  If matrices, they
% should be arranged as NCHANxLEN
%
% Assumes that spike times will not end up binned into same index for a
% given waveform (consistent with refractory period and small enough dt).
%

opts = inputParser();
opts.addParamValue('dt', 1);
opts.addParamValue('noisesigma', 0);
opts.addParamValue('nsamples', []);
opts.addParamValue('amps', {});
opts.addParamValue('shift', 0); % Useful for annoying case of spiketimes being 1-indexed...
opts.parse(varargin{:});
opts = opts.Results;

if isempty(opts.nsamples)
    opts.nsamples = ceil(max(cellfun(@(st) st(end), sts))/dt);
end


% Create the noise
maxchan = max(cellfun(@(w) size(w,2), waveforms));
if opts.noisesigma > 0
    trace = opts.noisesigma*randn(opts.nsamples, maxchan);
else
    trace = zeros(opts.nsamples, maxchan);
end


% Build up the output trace one waveform at a time
for i = 1:length(waveforms)
    st = sts{i} + opts.shift;
    thistrace = zeros(opts.nsamples, 1);
    spikeinds = round(st/opts.dt);
    spikeinds = spikeinds(spikeinds <= opts.nsamples);

    if isempty(opts.amps) || isempty(opts.amps{i})
        thistrace(spikeinds) = 1;
    else
        thistrace(spikeinds) = opts.amps{i};
    end
    
    nchan = size(waveforms{i}, 2);
    for c = 1:nchan
        trace(:,c) = trace(:,c) + conv(thistrace, waveforms{i}(:,c), 'same');
    end
end