function spiketimes = basic_sim_spikes(template, dx, trials, varargin)
% BASIC_SIM_SPIKES  Use spike rate template for simplistic spike simulation
%
% 2012-06 phli
%

opts = inputParser();
opts.addParamValue('baseline', 0);
opts.addParamValue('xstart', 1);
opts.addParamValue('xend', length(template));
opts.addParamValue('plot', nargout < 1);
opts.parse(varargin{:});
opts = opts.Results;

template = template + opts.baseline;
templatex = 1:length(template);
interpx = opts.xstart:dx:opts.xend;
interptemplate = interp1(templatex, template.*dx, interpx, 'linear', opts.baseline.*dx);

draws = rand(trials,length(interptemplate));
spikes = draws < repmat(interptemplate, [trials,1]);
[r c] = find(spikes);
spiketimes = interpx(c);


if opts.plot
    plot(spiketimes, r, 'k.');
    set(gca, 'XLim', [opts.xstart opts.xend]);
    clear spiketimes;
end