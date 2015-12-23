function varargout = ei_spike_usage(datarunorname)
% EI_SPIKE_USAGE    Calculate proportion of spikes used for EI calculation
% usage: ratio = ei_spike_usage(datarunorname)
%
% Input can be either datarun struct or string argument for LOAD_DATA.
% Output is the ratio of EI spikes to total spikes.  If no output is taken,
% then shows a histogram of these ratios.
%
% See also: LOAD_DATA
%
% 2011-09 phli
%


if ischar(datarunorname)
    datarun = load_data(datarunorname);
elseif isstruct(datarunorname)
    datarun = datarunorname;
else
    error('args');
end


if ~isfield(datarun, 'spikes')
    datarun = load_neurons(datarun);
end

if ~isfield(datarun, 'ei')
    datarun = load_ei(datarun, []);
end


spikes = cell2mat(collect(datarun.spikes, @length));
eispikes = datarun.ei.num_spikes';

ratio = eispikes ./ spikes;

if nargout > 0
    varargout{1} = ratio;
else
    hist(ratio, 0:0.1:1);
end
