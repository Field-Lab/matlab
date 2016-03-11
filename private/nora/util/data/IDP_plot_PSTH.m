function [spike_rate] = IDP_plot_PSTH(prepped_data, cell_to_plot, color)

if nargin == 2
    color = [0    0.4470    0.7410]; % default color 1
end
if length(color)==1
   default_colors = get(gca,'ColorOrder');
   color = default_colors(mod(color,7),:);
end

trials = size(prepped_data.testspikes, 1);
spikes = cell2mat(prepped_data.testspikes(:,cell_to_plot)); % spikes to frames
%duration = ceil(max(spikes));


bin = 1/119.5;
duration = 1200*bin;
bins = duration/bin;

spikes = floor(spikes/bin); % spikes to frames

spikes_per_bin = zeros(bins,1);
for i_bin = 1:bins
    spikes_per_bin(i_bin) = sum(spikes == i_bin);
end

spike_rate = conv(spikes_per_bin, gausswin(10), 'same');
spike_rate = spike_rate/(bin*trials*sum(gausswin(10)));
%spike_rate = spikes_per_bin/(bin*trials);

%time = (1:bins)*bin;
plot(spike_rate, 'Color' , color)

% set(gcf, 'Position', [ 500 500 800 200])

end