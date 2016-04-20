function [spike_rate, time] = IDP_plot_PSTH(prepped_data, cell_to_plot, varargin)

p = inputParser;
p.addParameter('color', [0    0.4470    0.7410]) % set color to 0 to just return the PSTH and not plot it
p.addParameter('bin', 1/119.5)
p.addParameter('smoothing', 1)
p.addParameter('duration', 10)
p.addParameter('plot_offset', 0)
p.parse(varargin{:});
color = p.Results.color;
bin = p.Results.bin;
smoothing = p.Results.smoothing;
duration = p.Results.duration;
plot_offset = p.Results.plot_offset;
clear p

if length(color)==1 && color ~= 0
   default_colors = get(gca,'ColorOrder');
   color = default_colors(mod(color,7)+1,:); 
end

trials = size(prepped_data.testspikes, 1);
spikes = cell2mat(prepped_data.testspikes(:,cell_to_plot));

%duration = ceil(max(spikes));
bins = floor(duration/bin);

% spike times to spikes per bin
spikes = floor(spikes/bin);
spikes_per_bin = zeros(bins,1);
for i_bin = 1:bins
    spikes_per_bin(i_bin) = sum(spikes == i_bin);
end

spike_rate = conv(spikes_per_bin, gausswin(smoothing), 'same'); % smooth the PSTH
spike_rate = spike_rate/(bin*trials*sum(gausswin(smoothing))); % spike count to spike rate

% plotting if that's your thing

time = (1:bins)*bin+plot_offset;
if length(color)==3
    plot(time, spike_rate, 'Color' , color)
    %set(gcf, 'Position', [ 500 500 800 200])
end

end