function proj_struct_new = spike_select_several_clusters(proj_struct,fig_to_use,num_clusters);

acf_range = 30; %in msec
cluster_colors = ['rgbmc']';
%cluster_colors = [0.2 0.2 0.2;0 1 0; 1 0.5 0; 0 0 1; 0 1 1;1 0 1; 1 .8 .8; 1 1 0];

% set up axes in the figure
proj_struct.plot_axes = spike_set_up_axes(fig_to_use,num_clusters);

%plot projections in projection_axes
proj_struct = spike_plot_density_gradient(proj_struct);

proj_struct.spike_projections = double(proj_struct.spike_projections);
for cc = 1:num_clusters

    % set up color for plot
    col = cluster_colors(mod(cc - 1,size(cluster_colors,2))+1,:);
    % have user select spikes
     s_i = spike_select_spikes(proj_struct,1,col);
     proj_struct.selected_indices{cc} = s_i{1};

    % note which spikes were selected
    selected_spike_waveforms = proj_struct.spike_waveforms(proj_struct.selected_indices{cc},:);
    selected_spike_times = proj_struct.spike_times(proj_struct.selected_indices{cc},1);

    % plot statistics about these spikes
    spike_compute_and_plot_acf(selected_spike_times,acf_range,proj_struct.plot_axes.acf_axis{cc});
    spike_plot_average_spike(selected_spike_waveforms,length(proj_struct.electrodes_for_projections),proj_struct.window_length,proj_struct.plot_axes.avg_spike_axis{cc},col);

end

proj_struct_new = spike_show_same_spikes_fancy(proj_struct,[],proj_struct,fig_to_use);
